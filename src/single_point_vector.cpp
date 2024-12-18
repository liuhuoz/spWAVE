#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

// Helper function: Calculate pairwise Euclidean distances between two matrices
Eigen::MatrixXd pairwise_dist(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2) {
    Eigen::MatrixXd mat1_sq = mat1.rowwise().squaredNorm(); // m x 1
    Eigen::MatrixXd mat2_sq = mat2.rowwise().squaredNorm(); // n x 1

    /*** 
    * Replicate and Extend the dimensions of squared norms
    * Eigen::MatrixXd requires same dimensions of two matrices to perform dot product. 
    ***/ 
    Eigen::MatrixXd mat1_sq_replicate = mat1_sq.replicate(1, mat2.rows()); // m x n
    Eigen::MatrixXd mat2_sq_replicate = mat2_sq.transpose().replicate(mat1.rows(), 1); // m x n

    // Compute cross-term (broadcasting)
    Eigen::MatrixXd cross_term = mat1 * mat2.transpose();

    // Final distance matrix
    Eigen::MatrixXd dist_mat = (mat1_sq_replicate + mat2_sq_replicate - 2 * cross_term).array().sqrt();
    return dist_mat;
}

// Helper function: transform Raw R matrix to Export Eigen MatrixXd
Eigen::MatrixXd mat_transform(
    const NumericMatrix& R_mat
){
    NumericMatrix mat_temp = Rcpp::as<NumericMatrix>(R_mat);
    Eigen::MatrixXd E_mat = Eigen::Map<Eigen::MatrixXd>(mat_temp.begin(), 
                                                        mat_temp.nrow(), 
                                                        mat_temp.ncol());
    return(E_mat);
}
//' Calculate single point vector through rcpp
//'
//' Calculate single molecular vector field from spatial coordinate and expression in one point through rcpp.
//'
//' @param X        Center point coordinates X.  
//' @param Y        Center point coordinates Y.  
//' @param spatial_expr     Spatial expression matrix, require x,y,gene expression.
//' @return         Dataframe contain Ex, Ey, U.
//' @seealso \code{\link{single_point_vector}}
//' @export 
// [[Rcpp::export]]
DataFrame single_point_vector_rcpp(
    const double& X,
    const double& Y,
    const NumericMatrix& spatial_expr
) {
    // Convert entire DataFrame (excluding rownames) to MatrixXd
    Eigen::MatrixXd full_mat = mat_transform(spatial_expr);

    // Separate coordinate columns and gene expression columns
    Eigen::MatrixXd coord_mat = full_mat.leftCols(2);     // First two columns: x, y
    Eigen::MatrixXd gene_expr_mat = full_mat.rightCols(full_mat.cols() - 2); // Remaining columns

    // Create source_mat (1x2 matrix with X and Y)
    Eigen::MatrixXd source_mat(1, 2);
    source_mat(0, 0) = X;
    source_mat(0, 1) = Y;

    // Calculate spatial distance matrix
    Eigen::MatrixXd spatial_dist_mat = pairwise_dist(coord_mat, source_mat);


    // Calculate spatial_sub_mat (vector differences)
    Eigen::MatrixXd spatial_sub_mat = source_mat.replicate(coord_mat.rows(), 1) - coord_mat;

    // Initialize result container
    Eigen::MatrixXd result(gene_expr_mat.cols(), 3); // Columns: Ex, Ey, U

    // Compute Ex, Ey, and U for each gene
    Eigen::VectorXd r_length_cubed = spatial_dist_mat.col(0).array().pow(3);

    for (int g = 0; g < gene_expr_mat.cols(); ++g) {
        Eigen::VectorXd gene_expr = gene_expr_mat.col(g);

        // Calculate Ex, Ey, U
        Eigen::VectorXd Ex = (spatial_sub_mat.col(0).array() / r_length_cubed.array()) * gene_expr.array();
        Eigen::VectorXd Ey = (spatial_sub_mat.col(1).array() / r_length_cubed.array()) * gene_expr.array();
        Eigen::VectorXd U = gene_expr.array() / spatial_dist_mat.col(0).array();

        // Handle NaN and Inf
        Ex = (Ex.array().isNaN() || Ex.array().isInf()).select(0, Ex);
        Ey = (Ey.array().isNaN() || Ey.array().isInf()).select(0, Ey);
        U = (U.array().isNaN() || U.array().isInf()).select(0, U);

        // Summarize results
        result(g, 0) = Ex.sum();
        result(g, 1) = Ey.sum();
        result(g, 2) = U.sum();
    }

    // Prepare output DataFrame
    // 提取列名为 CharacterVector
    CharacterVector col_genes = colnames(spatial_expr);

    // 将从第 3 列开始的基因名直接转换为 std::vector<std::string>
    std::vector<std::string> gene_names(
        col_genes.begin() + 2, // 从第 3 列开始（索引 2）
        col_genes.end()        // 到最后一列
    );

    Rcpp::DataFrame output = Rcpp::DataFrame::create(
        Named("gene") = gene_names,
        Named("Ex") = result.col(0),
        Named("Ey") = result.col(1),
        Named("U") = result.col(2)
    );

    return output;
}
