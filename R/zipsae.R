#' @title EBLUPs under Zero-Inflated Poisson Model
#' @description This function produces empirical best linier unbiased predictions (EBLUPs) for Zero-Inflated data and its Relative Standard Error. Small Area Estimation with Zero-Inflated Model (SAE-ZIP) is a model developed for Zero-Inflated data that can lead us to overdispersion situation. To handle this kind of situation, this model is created. The model in this package is based on Small Area Estimation with Zero-Inflated Poisson model proposed by Dian Christien Arisona (2018)<https://repository.ipb.ac.id/handle/123456789/92308>. For the data sample itself, we use combination method between Roberto Benavent and Domingo Morales (2015)<doi:10.1016/j.csda.2015.07.013> and Sabine Krieg, Harm Jan Boonstra and Marc Smeets (2016)<doi:10.1515/jos-2016-0051>.
#' @param data The data frame with vardir, response, and explanatory variables included with Zero-Inflated situation also.
#' @param vardir Sampling variances of direct estimations, if it is included in data frame so it is the vector with the name of sampling variances.if it is not, it is a data frame of sampling variance in order : \code{var1, cov12,.,cov1r,var2,cov23,.,cov2r,.,cov(r-1)(r),var(r)}
#' @param formula List of formula that describe the fitted model
#' @param PRECISION Limit of Fisher-scoring convergence tolerance. We set the default in \code{1e-4}
#' @param MAXITER Maximum number of iterations in Fisher-scoring algorithm. We set the default in \code{100}
#' @return This function returns a list of the following objects:
#'    \item{estimate}{A Vector with a list of EBLUP with Zero-Inflated Poisson model}
#'    \item{dispersion}{A list containing the following objects:}
#'      \itemize{
#'        \item rse : A dataframe with the values of relative square errors of estimation
#'      }
#'    \item{coefficient}{A list containing the following objects:}
#'      \itemize{
#'        \item lambda : The estimator of model based on Non-Zero data
#'        \item omega : The estimator of model based Complete Data
#'      }
#' @examples
#' ##load the dataset in package
#' data(dataSAEZIP)
#'
#' ##Extract the vardir (sampling error)
#' dataSAEZIP$vardir -> sError
#'
#' ##Compute the data with SAE ZIP model
#' formula = (y~x1)
#' zipsae(data = dataSAEZIP, vardir = sError, formula) -> saezip
#'
#' saezip$estimate        #to see the result of Small Area Estimation with Zero-Inflated Model
#' saezip$dispersion$rse  #to see the relative standard error from the estimation
#' saezip$coefficient$lambda   #to see the estimator which is gained from the non-zero compilation data
#' saezip$coefficient$omega   #to see the estimator which is gained from the complete compilation data.
#'
#' head(saezip)
#'
#' @export zipsae
#' @import stats
zipsae <- function(data, vardir, formula, PRECISION = 1e-04, MAXITER = 100 ){
  if(missing(data)){
    stop("Data argument is missing")
  }
  if (length(formula)==0){
    stop("Please insert your formula for fitted model")
  }
  if (missing(vardir)){
    stop("vardir argument is missing")
  }
  hasil <- list(estimate = NA,
                dispersion = list(rse = NA),
                coefficient = list(lambda = NA,
                                 omega = NA
                                 )
            )
  reml <- function(X, Y, vardir, MAXITER){
    PRECISION = 1e-04
    REML_tmp <- 0
    REML_tmp[1] <- mean(vardir)
    k <- 0
    acc <- PRECISION + 1
    while ((acc > PRECISION) & (k < MAXITER)) {
      k <- k + 1
      Psi <- 1/(REML_tmp[k] + vardir)
      PsiXt <- t(Psi * X)
      Q <- solve(PsiXt %*% X)
      P <- diag(Psi) - t(PsiXt) %*% Q %*% PsiXt
      Py <- P %*% Y
      s <- (-0.5) * sum(diag(P)) + 0.5 * (t(Py) %*% Py)
      Fis <- 0.5 * sum(diag(P %*% P))
      REML_tmp[k + 1] <- REML_tmp[k] + s/Fis
      acc <- abs((REML_tmp[k + 1] - REML_tmp[k])/REML_tmp[k])
    }
    REML <- max(REML_tmp[k + 1], 0)
    return(REML)
  }
  m_data = model.frame(formula, na.action = na.omit ,data)
  Xz <- model.matrix(formula, m_data)
  Xz <- as.numeric(Xz[,-1])
  Xnz <- Xz
  Y <- m_data[, 1]
  ToOmit <- c()
  dataLength <- nrow(data)
  for (i in 1:dataLength) {
    if (Y[i] == 0) {
      ToOmit <- c(ToOmit, i)
    }
  }
  if (length(ToOmit) != 0) {
    Xnz <- Xz[-ToOmit]
    Ynz <- Y[-ToOmit]
    vardirNz <- vardir[-ToOmit]
  }else{
    Xnz <- Xz
    Ynz <- Y
    vardirNz <- vardir
  }
  b.REMLz  <- reml(Xz, Y, vardir, MAXITER)
  a.REMLnz <- reml(Xnz, Ynz, vardirNz, MAXITER)
  Vnz <- 1/(a.REMLnz + vardirNz)
  Vz <- 1/(b.REMLz + vardir)
  Xnz_tmp <- t(Vnz*Xnz)
  Xz_tmp <- t(Vz*Xz)
  Qnz <- solve(Xnz_tmp%*%Xnz)
  Qz <- solve(Xz_tmp%*%Xz)
  beta.nz <- Qnz%*%Xnz_tmp%*%Ynz
  beta.z <- Qz%*%Xz_tmp%*%Y
  X_beta.nz <- Xnz%*%beta.nz
  X_beta.z <- Xz%*%beta.z
  residNz <- Ynz - X_beta.nz
  residZ <- Y - X_beta.z
  lamda = exp((Xz%*%beta.nz) + as.vector(a.REMLnz*(Vnz%*%residNz)))
  logit = exp(X_beta.z + as.vector(b.REMLz*(Vz%*%residZ)))
  omega = (logit/(1 + logit))
  delta = 1-omega
  est <- lamda*delta
  est - mean(est) -> diff_tmp
  (diff_tmp)^2 -> diff_tmp_squared
  sqrt(diff_tmp_squared)/est -> rse_tmp
  hasil$estimate <- est
  hasil$dispersion$rse <- rse_tmp
  hasil$coefficient$lambda <- lamda
  hasil$coefficient$omega <- omega
  return(hasil)
}
