#**********************************

#**************************
#** Tensor Decomposition **
#**************************


#' prepare tensor data for CP decomposition
#'
#' This function pre-processes the LR pair field data from a spWAVE object
#'
#' @param spwave spWAVE object containing LR pair field data
#' @param low_quantile quantile threshold to filter low expression values. Default is 0.1 (10%).
#'
#' @return log-normalized filtered tensor array for CP decomposition
#' @export
prep_TD_data <- function(spwave,low_quantile=0.1){
  # create tensor data
  tensor_data <- spwave@LR_pair_field

  tensor_data_log <- array(data=NA,
    dim=c(dim(tensor_data)[1],3,dim(tensor_data)[3]),
    dimnames = list(dimnames(tensor_data)[[1]],c("logNorm","logEx","logEy"),dimnames(tensor_data)[[3]])
  )

  tensor_data_log[,"logNorm",] <-
    log(tensor_data[,"Ex",]^2+tensor_data[,"Ey",]^2)/2 #log 歪除2等于取了开方，节省计算开销
  tensor_data_log[,"logEx",] <-
    abs(tensor_data_log[,"logNorm",])*cos(atan2(tensor_data[,"Ey",],tensor_data[,"Ex",]))
  tensor_data_log[,"logEy",] <-
    abs(tensor_data_log[,"logNorm",])*sin(atan2(tensor_data[,"Ey",],tensor_data[,"Ex",]))

  # filter low expression
  tensor_norm <-
    apply(tensor_data,3,function(x){
      sqrt(x[,"Ex"]^2+x[,"Ey"]^2)
    })

  low_cutoff <- tensor_norm %>% as.vector %>% quantile(probs=c(low_quantile))

  cond <- tensor_data_log[, "logNorm", ] < log(low_cutoff)
  cond_expanded <- array(cond, dim = c(dim(cond), ncol(tensor_data_log)))
  cond_expanded <- aperm(cond_expanded, perm = c(1, 3, 2))
  tensor_data_log[cond_expanded] <- 0
  return(tensor_data_log[,2:3,])
}


#' Repeated Tensor CP Decomposition
#'
#' This function performs CP decomposition on a tensor multiple times for a range of ranks
#'
#' @param tensor Input tensor array for CP decomposition
#' @param num_rep Number of repetitions for each rank. Default is 10.
#' @param rank_range Range of ranks to test for CP decomposition. Default is 10:15.
#' @param random_seed Seed for random number generation to ensure reproducibility. Default is 42.
#' @param const character vextor, Constraints for each mode during CP decomposition. Default is c("uncons","uncons","orthog").
#' see multiway::parafac and CMLS::cmls for more details and available options.
#' @param corcondia_cutoff Minimum corcondia value to accept a decomposition. Default is 10.
#' @param parallel Logical indicating whether to use parallel processing. Default is TRUE.
#' @param cl Cluster object for parallel processing. Default is NULL.
#'
#' @importFrom multiway parafac
#' @importFrom multiway corcondia
#'
#' @seealso \code{\link[multiway]{parafac}}
#' @seealso \code{\link[CMLS]{cmls}}
#' @seealso \code{\link[multiway]{corcondia}}
#'
#' @return nested list of CP decomposition results for each rank
#' @export
repeat_tensor_CPD <- function(
  tensor,
  num_rep=10,
  rank_range=10:15,
  random_seed=42,
  const=c("uncons","uncons","orthog"),
  corcondia_cutoff=10,
  parallel=TRUE,cl=NULL){
  set.seed(random_seed)

  if (!requireNamespace("multiway", quietly = TRUE)) {
    stop(
      "Package \"multiway\" must be installed to use this function.",
      call. = FALSE
    )
  }


  if(parallel){
    if (!requireNamespace("parallel", quietly = TRUE)) {
    stop(
      "Package \"parallel\" must be installed to use this function.",
      call. = FALSE
    )
  }
    if(!is.null(cl)){
      clusterEvalQ(cl,library(multiway))
    }else{
      stop("Please provide a valid cluster object 'cl' for parallel processing.")
    }
  }

  CP_decomp <- list()
  loop_idx <- seq_along(rank_range)
  for(i in loop_idx){
    R=rank_range[i]

    cat("Processing rank", R, ":\n")

    rep_results <- list()
    for(j in 1:num_rep) {
      cat("Repeat", j, "...\n")
      while(length(rep_results) < j){
        temp_result <- multiway::parafac(tensor,nfac=R,
            const=const,
            parallel=parallel,cl=cl)
        corcondia_temp <- multiway::corcondia(tensor, temp_result)
        print(corcondia_temp)
        if(corcondia_temp > corcondia_cutoff){
          rep_results[[j]] <- temp_result
        }else{
          cat("corcondia below", corcondia_cutoff, "try again...\n")
        }
      }
    }
    CP_decomp[[i]] <- rep_results
  }
  names(CP_decomp) <- paste0("Rank_", rank_range)
  return(CP_decomp)
}
