#' ace_sampl_tip
#' 
#' ace function with sampling
#' Amended version of the ape::ace fucntion incorporating sampling rate estimates at the tips when calculating the movement matrix between demes.
#' Set up to work for a discrete charater (location) with the matrix exponentiation technique. Eigenvalue and continuous data methods have been removed from the function. 
#' The default model has been set to ARD, all rates different, as this adjustment has been tested on a non-symetric matrix only to provide a useful estimate of migration between locations. 
#' 
#' 
#' @param x a named vector or factor providing the location of the tips included in the input tree. A factor is recommended for inference and sampling rate incorporation as the order of the locations can be set and unobserved locations can be included.
#' @param phy input phylogeny over which to estimate rates. Object of class "phylo".
#' @param sampling a vector of sampling rates by location. These will be converted to relative rates in the function.
#' @param CI a logical specifying whether to return the likelihood of the different discrete states at the root.
#' @param model character specifying the matrix for the model. ARD for all rates different (eg. "\code{matrix(c(0, 1, 2, 3, 0, 4, 5, 6, 0), 3)}"), SYM for symetrical (eg. "\code{matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3)}") and ER for equal rates (eg. "\code{matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3)}").
#' @param kappa a positive value giving the exponent transformation of the branch lengths. Similar to the correction but same across all branches.
#' @param ip initial value(s) for the ML estimation.
#' @param use.expm logical specifying whether to to use the expm package to compute the matrix exponential (only tested on this). Alternative is equivalent ape matexpo function.
#' @param marginal set true for marginal rather than joint reconstruction of ancestral states (not tested). This marginal method is the traditional reconstruction method where only descendents are used to estimate the likelihood at a node.
#' @param lwbd lower bound for the ML estimation.
#' @param upbd upper bound for the ML estimation.
#'
#' @return object of class 'ace'
#' @export
#'
#' @examples
ace_sampl_tip <- function (x
                           , phy
                           , sampling
                           , CI = FALSE
                           , model = "ARD_rel"
                           , kappa = 1
                           , ip = c(3,0) # can specify a vector with different starting values for a and d params
                           , use.expm = TRUE
                           , marginal = FALSE
                           , lwbd = -1e+50       ## added explicit upper and lower bound parameters, can't go to 0 with relative rates as in the original code
                           , upbd = 1e+50) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("tree has no branch lengths")
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("\"phy\" is not rooted AND fully dichotomous.")
  if (length(x) != nb.tip) 
    stop("length of phenotypic and of phylogenetic data do not match.")
  if (!is.null(names(x))) {
    if (all(names(x) %in% phy$tip.label)) 
      x <- x[phy$tip.label]
    else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
  }
  obj <- list()
  if (kappa != 1) 
    phy$edge.length <- phy$edge.length^kappa
  
  if (any(phy$edge.length <= 0))
    stop("some branches have length zero or negative")
  if (!is.factor(x)) 
    x <- factor(x)
  nl <- nlevels(x)
  lvls <- levels(x)
  x <- as.integer(x)
  if (is.character(model)) {
    rate <- matrix(NA, nl, nl)
    switch(model, ER = np <- rate[] <- 1, ARD_rel = {
      np <- nl * (nl - 1)
      sel <- col(rate) < row(rate) ## adapted from the SYM model, additional param for relative added later
      rate[sel] <- 1:(np/2) 
      rate <- t(rate)
      rate[sel] <- 1:(np/2)
      },
      SYM = {
      np <- nl * (nl - 1)/2
      sel <- col(rate) < row(rate)
      rate[sel] <- 1:np
      rate <- t(rate)
      rate[sel] <- 1:np
    })
  }
  else {
    if (ncol(model) != nrow(model)) 
      stop("the matrix given as 'model' is not square")
    if (ncol(model) != nl) 
      stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
    rate <- model
    np <- max(rate)
  }
  index.matrix <- rate
  tmp <- cbind(1:nl, 1:nl)
  index.matrix[tmp] <- NA
  rate[tmp] <- 0
  rate[rate == 0] <- np/2+1 ## needs to be one bigger than the maximum index
  liks <- matrix(0, nb.tip + nb.node, nl)
  TIPS <- 1:nb.tip
  liks[cbind(TIPS, x)] <- 1
  if (anyNA(x)) 
    liks[which(is.na(x)), ] <- 1
  phy <- reorder(phy, "postorder")
  Q <- matrix(0, nl, nl)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  if (!requireNamespace("expm", quietly = TRUE) && 
      use.expm) {
    warning("package 'expm' not available; using function 'matexpo' from 'ape'")
    use.expm <- FALSE
  }
  E <- if (use.expm) 
    expm::expm
  else ape::matexpo
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  dev <- function(p, output.liks = FALSE) {
    if (is.null(sampling)){
      if (any(is.nan(p)) || any(is.infinite(p))) 
        return(1e+50)
      comp <- numeric(nb.tip + nb.node)
      p <- exp(p)
      
      #### Relative method ####
      a <- p[1:(length(p)/2)]
      d <- p[((length(p)/2)+1):length(p)]
      Q[] <- c(a, 0)[rate]
      Q[upper.tri(Q)] <- a*d

      diag(Q) <- -rowSums(Q)
      
      decompo <- eigen(Q)
      lambda <- decompo$values
      GAMMA <- decompo$vectors
      invGAMMA <- solve(GAMMA)
      for (i in seq(from = 1, by = 2, length.out = nb.node)) {
        j <- i + 1L
        anc <- e1[i]
        des1 <- e2[i]
        des2 <- e2[j]
        v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% 
          invGAMMA %*% liks[des1, ]
        v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% 
          invGAMMA %*% liks[des2, ]
        v <- v.l * v.r
        comp[anc] <- sum(v)
        liks[anc, ] <- v/comp[anc]
      }
      if (output.liks) 
        return(liks[-TIPS, ])
        dev <- -2*sum(log(comp[-TIPS]))

      if (is.na(dev)) 
        Inf
      
      else dev
    }else{
      if (any(is.nan(p)) || any(is.infinite(p))) 
        return(1e+50)
      comp <- numeric(nb.tip + nb.node)
        p <- exp(p) 
      
      #### Relative method ####
      a <- p[1:(length(p)/2)]
      d <- p[((length(p)/2)+1):length(p)]
      #d <- exp(d)
      Q[] <- c(a, 0)[rate]
      Q[upper.tri(Q)] <- a*d
      
      #### Changes in correction ####
      # These two optins are tested in the text of the thesis, sampling/sum(sampling) was found to perform best
      probs <- sampling
      probs <- (sampling/sum(sampling))
      
      
      
      
      Q2 <- Q%*%diag(probs)
      diag(Q) <- -rowSums(Q) 
      diag(Q2)<- -rowSums(Q2)
      decompo <- eigen(Q)
      decompo2 <- eigen(Q2)
      lambda <- decompo$values
      lambda2 <- decompo2$values
      GAMMA <- decompo$vectors
      GAMMA2 <- decompo2$vectors
      invGAMMA <- solve(GAMMA)
      invGAMMA2 <- solve(GAMMA2)
      
      for (i in seq(from = 1, by = 2, length.out = nb.node)) {
        j <- i + 1L
        anc <- e1[i]
        des1 <- e2[i]
        des2 <- e2[j]
        if(des1%in%TIPS){
           v.l <- GAMMA2 %*% diag(exp(lambda2 * EL[i])) %*% 
             invGAMMA2 %*% liks[des1, ]
         }else{
         v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% 
            invGAMMA %*% liks[des1, ]
         }
         if(des2%in%TIPS){
           v.r <- GAMMA2 %*% diag(exp(lambda2 * EL[j])) %*% 
             invGAMMA2 %*% liks[des2, ]
         }else{
          v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% 
             invGAMMA %*% liks[des2, ]
         }
        v <- v.l * v.r
        comp[anc] <- sum(v)
        liks[anc, ] <- v/comp[anc]
      }
      if (output.liks) 
        return(liks[-TIPS, ])
        dev <- -2*sum(log(comp[-TIPS]))
      if (is.na(dev)) 
        Inf
      
      
      else dev
    }
  }
  
  
out <- try(optim(par = rep(ip, length.out = np), fn = function(p) dev(p), method = c('L-BFGS-B'), lower = c(-7,-10), upper = c(7,10), hessian = T, control = c(maxit = 10000)))

  if(class(out) == "try-error"){
    warning("convergence failed")
    obj$loglik <- NaN  
    obj$rates <- rep(NaN, np)
    obj$se <- rep(NaN, np)
  } else {
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    
      obj$loglik <- -out$value/2
      obj$rates <- out$par
     obj$se <- ape:::.getSEs(out)
  }
  obj$index.matrix <- index.matrix
  
  if (CI) {
    lik.anc <- dev(obj$rates, TRUE)
    if (!marginal) {
      Q[] <- c(obj$rates, 0)[rate]
      diag(Q) <- -rowSums(Q)
      for (i in seq(to = 1, by = -2, length.out = nb.node)) {
        anc <- e1[i] - nb.tip
        des1 <- e2[i] - nb.tip
        if (des1 > 0) {
          P <- matexpo(Q * EL[i])
          tmp <- lik.anc[anc, ]/(lik.anc[des1, ] %*% P)
          lik.anc[des1, ] <- (tmp %*% P) * lik.anc[des1, ]
        }
        j <- i + 1L
        des2 <- e2[j] - nb.tip
        if (des2 > 0) {
          P <- matexpo(Q * EL[j])
          tmp <- lik.anc[anc, ]/(lik.anc[des2, ] %*%  P)
          lik.anc[des2, ] <- (tmp %*% P) * lik.anc[des2, ]
        }
        lik.anc <- lik.anc/rowSums(lik.anc)
      }
    }
    colnames(lik.anc) <- lvls
    obj$lik.anc <- lik.anc
  }
  obj$call <- match.call()
  class(obj) <- "relace"
  obj
}






## changes to print S3 classes

logLik.relace <- function(object, ...) object$loglik

deviance.relace <- function(object, ...) -2*object$loglik

AIC.relace <- function(object, ..., k = 2)
{
  if (is.null(object$loglik)) return(NULL)
  ## Trivial test of "type"; may need to be improved
  ## if other models are included in ace(type = "c")
  np <- if (!is.null(object$sigma2)) 1 else length(object$rates)
  -2*object$loglik + np*k
}

### by BB:
anova.relace <- function(object, ...)
{
  X <- c(list(object), list(...))
  df <- lengths(lapply(X, "[[", "rates"))
  ll <- sapply(X, "[[", "loglik")
  ## check if models are in correct order
  dev <- c(NA, 2*diff(ll))
  ddf <- c(NA, diff(df))
  table <- data.frame(ll, df, ddf, dev,
                      pchisq(dev, ddf, lower.tail = FALSE))
  dimnames(table) <- list(1:length(X), c("Log lik.", "Df",
                                         "Df change", "Resid. Dev",
                                         "Pr(>|Chi|)"))
  structure(table, heading = "Likelihood Ratio Test Table",
            class = c("anova", "data.frame"))
}

print.relace <- function(x, digits = 4, ...)
{
  cat("\n    Relative Ancestral Character Estimation\n\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  if (!is.null(x$loglik))
    cat("    Log-likelihood:", x$loglik, "\n\n")
  if (!is.null(x$resloglik))
    cat("    Residual log-likelihood:", x$resloglik, "\n\n")
  ratemat <- x$index.matrix
  if (is.null(ratemat)) { # to be improved
    class(x) <- NULL
    x$resloglik <- x$loglik <- x$call <- NULL
    print(x)
  } else {
    dimnames(ratemat)[1:2] <- dimnames(x$lik.anc)[2]
    cat("Rate index matrix:\n")
    print(ratemat, na.print = ".")
    cat("\n")
    npar <- length(x$rates)
    estim <- data.frame((1:npar), round(x$rates, digits), round(x$se, digits))
    cat("Parameter estimates:\n")
    names(estim) <- c("rate index", "estimate", "std-err")
    print(estim, row.names = FALSE)
    if (!is.null(x$lik.anc)) {
      cat("\nScaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):\n")
      print(x$lik.anc[1, ])
    }
  }
}