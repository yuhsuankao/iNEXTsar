H_hat_inc_wor = function(y, T_star){

  nT = y[1]
  y = y[-1]
  y = y[y > 0]
  U = sum(y)

  rho = nT/T_star
  Q1 = sum(y == 1)
  Q2 = sum(y == 2)

  js = sort(unique(y)) #To only calculate unique yi once in the combination
  Qj = table(y)
  # k = 1:(nT-1)
  # term1 = sapply(js, function(x) x/nT * exp(lchoose(nT-x,k)-lchoose(nT-1,k)) ) %>%
  #   apply(1, '*', as.vector(Qj)) %>% t %>% apply(2, '*', 1/k*(T_star-k)/T_star) %>% sum

  temp = sapply(js, function(z){

    k = 1:(nT-z);k = k[k>0]
    sum(1/k*(T_star-k)/T_star*z/nT*exp(lchoose(nT-z,k)-lchoose(nT-1,k)))

  })
  term1 = sum(temp * Qj)

  B_star = if(Q1==0&Q2==0){ 1 } else {((1-rho)*2*Q2 + rho*Q1) / ((nT-1)*Q1 + 2*Q2)}

  if (B_star == 1) {term2 = 0} else {
    r = 1:(nT-1)
    term2 =  (T_star-nT)/T_star * (Q1/nT) * (1-B_star)^(-nT+1) * (-log(B_star) - sum(1/r * (1-B_star)^r))
  }

  H0_hat = term1 + term2
  H_hat = nT/U * H0_hat + log(U/nT)
  return( H_hat )

}
S_hat_inc_wor = function(y, T_star) {

  nT = y[1]
  y = y[-1]; y = y[y > 0]

  rho = nT/T_star
  Q1 = sum(y == 1)
  Q2 = sum(y == 2)

  Sorb = length(y)
  Q0_hat = if(Q1==0&Q2==0) { 0 } else { Q1^2 / ( nT/(nT-1)*2*Q2 + rho/(1-rho)*Q1 ) }
  return(Sorb+Q0_hat)

}
Dqhat.Sam_inc = function (y, T_star, q, t) {
  nT = y[1]
  y = y[-1]; y = y[y > 0]
  U = sum(y)
  rho = nT/T_star
  Q1 = sum(y == 1)
  Q2 = sum(y == 2)

  Q0_hat = if(Q1==0&Q2==0 | rho==1 | nT==1) { 0 } else { Q1^2 / ( nT/(nT-1)*2*Q2 + rho/(1-rho)*Q1 ) }
  M0_bar = if(Q0_hat==0 | rho==0) { 0 } else {(1/rho-1)*Q1/Q0_hat}
  M1_bar = if(Q1==0&Q2==0 | nT==1) { 0 } else if (Q1==0){ T_star-nT+2 } else {1 + (T_star-nT)*2*Q2 / ((nT-1)*Q1)}
  sum_pi_square = if(nT==1 | nT==0) { 0 } else {(sum(y*(y-1)) + U*rho) / (nT*(nT-1) + nT*rho)}
  sum_pi = U/nT

  Qk.hat = function(y, nT, t) {
    y <- y[y > 0]
    Qj <- table(y)
    js <- as.numeric(names(Qj))

    Sub <- function(k) {
      Qj_k <- Qj[js >= k]
      js_k <- js[js >= k]
      ifelse(length(js_k) == 0, 0, sum(exp(lchoose(js_k,k) + lchoose(nT - js_k, t - k) - lchoose(nT,t)) * Qj_k))
    }
    sapply(1:t, Sub)
  }

  delta0_hat = function(y, nT, t) {
    Sub <- function(t) {
      if (t <= nT) {
        Fun <- function(y) {
          if (y <= (nT - t))
            exp(lgamma(nT - y + 1) + lgamma(nT - t +1) - lgamma(nT - y - t + 1) - lgamma(nT +1))
          else 0
        }
        sum(1 - sapply(y, Fun))
      }
      else {
        Sobs = sum(y > 0)

        delta0_hat = ifelse(rho==1, Sobs, Sobs + Q0_hat*(1-(1-((t-nT)/(T_star-nT)))^M0_bar))
        return(delta0_hat)
      }
    }
    sapply(t, Sub)
  }

  delta1_hat = function(y, nT, t) {

    Sub <- function(t) {
      if (t <= nT) {
        k <- 1:t
        Ut.hat <- t/nT * U
        exp(-sum(k/Ut.hat * log(k/Ut.hat) * Qk.hat(y,nT, t)))
      }
      else {
        if(rho==1){exp(-sum(y/U*log(y/U)))} else{
          ( exp(H_hat_inc_wor(c(nT,y), T_star)) - exp(-sum(y/U*log(y/U))) ) * (1-(1-(t-nT)/(T_star-nT))^M0_bar) + exp( -sum(y/U*log(y/U)) )
        }}
    }
    sapply(t, Sub)
  }

  delta2_hat = function(y, nT, t) {

    Sub <- function(t){

      if(t<=nT) {1/(1/t * nT/U + (1 - 1/t) * sum(y *(y - 1)/U^2/(1 - 1/nT)))

      }
      else{ ( (t*U/nT)^(-2) * ((1-t/T_star)*t*sum_pi - (1-t/T_star)*t*sum_pi_square + t^2*sum_pi_square) )^(-1) }

    }
    sapply(t, Sub)
  }

  if (q == 0)
    delta0_hat(y, nT, t)
  else if (q == 1)
    delta1_hat(y, nT, t)
  else if (q == 2)
    delta2_hat(y, nT, t)
}
Chat.Sam_inc = function(y, T_star, t){

  nT = y[1]
  y = y[-1]; y = y[y > 0]
  U = sum(y)
  Q1 = sum(y == 1)
  Q2 = sum(y == 2)

  M1_bar = if(Q1==0&Q2==0 | nT==1) { 0 } else if (Q1==0){ T_star-nT+2 } else {1 + (T_star-nT)*2*Q2 / ((nT-1)*Q1)}

  Sub = function(t){

    if (t<nT){

      # uni_yi = sort(unique(y)) #To only calculate unique yi once in the combination
      # Qi = as.vector(table(y))
      #
      # 1 - ((T_star-t)/T_star)*sum( uni_yi/U * exp( lchoose(nT-uni_yi,t)-lchoose(nT-1,t) ) * Qi )

      1-(T_star-t)/T_star*sum(exp(lchoose(nT-y,t)-lchoose(nT-1,t))*y/U)

    } else {

      if(nT==T_star) { 1 } else { 1 - ((T_star-nT)/T_star)*Q1/U*(1-(t-nT)/(T_star-nT))^M1_bar }

    }
  }
  sapply(t, Sub)
}

iNEXT.Sam_sar_inc = function (Spec, t = NULL, q = 0, endpoint = min(2*max(Spec), T_star), knots = 40, se = TRUE, nboot = 200, conf = 0.95, T_star) {
  if (which.max(Spec) != 1)
    stop("invalid data structure!, first element should be number of sampling units")
  nT <- Spec[1]
  if (is.null(t)) {
    if (endpoint <= nT) {
      t <- floor(seq(1, endpoint, length.out = floor(knots)))
    }
    else {
      t <- c(floor(seq(1, nT - 1, length.out = floor(knots/2) - 1)), nT, floor(seq(nT + 1, to = endpoint, length.out = floor(knots/2))))
    }
    t <- c(1, t[-1])
  }
  else if (is.null(t) == FALSE) {
    if (max(t) > nT & length(t[t == nT]) == 0)
      t <- c(t, nT - 1, nT, nT + 1)
    t <- sort(t)
  }
  Dq.hat <- Dqhat.Sam_inc(Spec, T_star, q, t)
  C.hat <- Chat.Sam_inc(Spec, T_star, t)
  if (se == TRUE & nboot > 0 & length(Spec) > 2) {
    Prob.hat <- EstiBootComm.Sam(Spec)
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot,
                                                      nT, p)))
    Abun.Mat <- matrix(c(rbind(nT, Abun.Mat)), ncol = nboot)
    tmp <- which(colSums(Abun.Mat) == nT)
    if (length(tmp) > 0)
      Abun.Mat <- Abun.Mat[, -tmp]
    if (ncol(Abun.Mat) == 0) {
      out <- cbind(t = t, qD = Dq.hat, SC = C.hat)
      warning("Insufficient data to compute bootstrap s.e.")
    }
    else {
      error <- qnorm(1 - (1 - conf)/2) * apply(apply(Abun.Mat,
                                                     2, function(y) Dqhat.Sam(y, q, t)), 1, sd, na.rm = TRUE)
      left <- Dq.hat - error
      right <- Dq.hat + error
      left[left <= 0] <- 0
      error.C <- qnorm(1 - (1 - conf)/2) * apply(apply(Abun.Mat,
                                                       2, function(y) Chat.Sam(y, t)), 1, sd, na.rm = TRUE)
      left.C <- C.hat - error.C
      right.C <- C.hat + error.C
      left.C[left.C <= 0] <- 0
      right.C[right.C >= 1] <- 1
      out <- cbind(t = t, qD = Dq.hat, qD.LCL = left,
                   qD.UCL = right, SC = C.hat, SC.LCL = left.C,
                   SC.UCL = right.C)
    }
  }
  else {
    out <- cbind(t = t, qD = Dq.hat, SC = C.hat)
  }
  out <- data.frame(out)
  out$method <- ifelse(out$t < nT, "interpolated", ifelse(out$t ==
                                                            nT, "observed", "extrapolated"))
  out$order <- q
  id <- match(c("t", "method", "order", "qD", "qD.LCL", "qD.UCL",
                "SC", "SC.LCL", "SC.UCL"), names(out), nomatch = 0)
  out <- out[, id]
  return(out)
}

iNEXT_sar_inc = function (x, q = 0, datatype = "abundance", size = NULL, endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = 50, T_star) {

  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence") {
    stop("datatype=\"incidence\" was no longer supported after v2.0.8, \n         please try datatype=\"incidence_freq\".")
  }
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  Fun <- function(x, q) {
    x <- as.numeric(unlist(x))
    if (datatype == "abundance") {
      if (sum(x) == 0)
        stop("Zero abundance counts in one or more sample sites")
      out <- iNEXT.Ind(Spec = x, q = q, m = size, endpoint = ifelse(is.null(endpoint),
                                                                    2 * sum(x), endpoint), knots = knots, se = se,
                       nboot = nboot, conf = conf)
    }
    if (datatype == "incidence") {
      t <- x[1]
      y <- x[-1]
      if (t > sum(y)) {
        warning("Insufficient data to provide reliable estimators and associated s.e.")
      }
      if (sum(x) == 0)
        stop("Zero incidence frequencies in one or more sample sites")
      out <- iNEXT.Sam_sar_inc(Spec = x, q = q, t = size, endpoint = ifelse(is.null(endpoint), min(2*max(x), T_star), endpoint), knots = knots, se = se,
                               nboot = nboot, conf = conf, T_star)
    }
    out
  }
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a postive value/vector of numeric object")
  if (min(q) < 0) {
    warning("ambigous of order q, we only compute postive q")
    q <- q[q >= 0]
  }
  if (class(x) == "numeric" | class(x) == "integer" | class(x) ==
      "double") {
    out <- do.call("rbind", lapply(q, function(q) Fun(x,
                                                      q)))
    out[, -(1:3)] <- round(out[, -(1:3)], 3)
    # index <- rbind(as.matrix(ChaoSpecies_wor_inc(x, datatype, conf, T_star=T_star)),
    #                as.matrix(ChaoEntropy_wor_inc(x, datatype, transform = TRUE,
    #                                              conf, T_star=T_star)), as.matrix(EstSimpson_wor_inc(x, datatype, transform = TRUE,
    #                                                                                                  conf, T_star=T_star)))
    # rownames(index) <- c("Species Richness", "Shannon diversity",
    #                      "Simpson diversity")
  }
  else if (class(x) == "matrix" | class(x) == "data.frame") {
    out <- apply(as.matrix(x), 2, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,
                                                        q)))
      tmp[, -(1:3)] <- round(tmp[, -(1:3)], 3)
      tmp
    })
    # arr <- array(0, dim = c(3, 5, ncol(x)))
    # arr[1, , ] <- t(as.matrix(ChaoSpecies_wor_inc(x, datatype, conf, T_star=T_star)))
    # arr[2, , ] <- t(as.matrix(ChaoEntropy_wor_inc(x, datatype, transform = TRUE,
    #                                               conf, T_star=T_star)))
    # arr[3, , ] <- t(as.matrix(EstSimpson_wor_inc(x, datatype, transform = TRUE,
    #                                              conf, T_star=T_star)))
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity",
    #                         "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.",
    #                         "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3, 1))
    # index <- dcast(as.data.frame(index), formula = Var1 +
    #                  Var2 ~ Var3, value.var = "Freq")
    # colnames(index) <- c("Site", "Diversity", "Observed",
    #                      "Estimator", "s.e.", "LCL", "UCL")
    # #just for no bootstrap
    # index = index[,1:4]
  }
  else if (class(x) == "list") {
    out <- lapply(x, function(x) {
      tmp <- do.call("rbind", lapply(q, function(q) Fun(x,
                                                        q)))
      tmp[, -(1:3)] <- round(tmp[, -(1:3)], 3)
      tmp
    })
    # arr <- array(0, dim = c(3, 5, length(x)))
    # arr[1, , ] <- t(as.matrix(ChaoSpecies_wor_inc(x, datatype, conf, T_star=T_star)))
    # arr[2, , ] <- t(as.matrix(ChaoEntropy_wor_inc(x, datatype, transform = TRUE,
    #                                               conf, T_star=T_star)))
    # arr[3, , ] <- t(as.matrix(EstSimpson_wor_inc(x, datatype, transform = TRUE,
    #                                              conf, T_star=T_star)))
    # dimnames(arr)[[3]] <- names(x)
    # dimnames(arr)[[1]] <- c("Species richness", "Shannon diversity",
    #                         "Simpson diversity")
    # dimnames(arr)[[2]] <- c("Observed", "Estimator", "Est_s.e.",
    #                         "Lower_CI", "Upper_CI")
    # index <- ftable(arr, row.vars = c(3, 1))
    # index <- dcast(as.data.frame(index), formula = Var1 +
    #                  Var2 ~ Var3, value.var = "Freq")
    # colnames(index) <- c("Site", "Diversity", "Observed",
    #                      "Estimator", "s.e.", "LCL", "UCL")
    # #just for no bootstrap
    # index = index[,1:4]
  }
  else {
    stop("invalid class of x, x should be a object of numeric, matrix, data.frame, or list")
  }
  info <- DataInfo(x, datatype)
  # z <- list(DataInfo = info, iNextEst = out, AsyEst = index)
  z <- list(DataInfo = info, iNextEst = out)
  class(z) <- c("iNEXT")
  return(z)
}



BootstrapFun = function (x, FunName, datatype, B) {
  if (!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
    stop("invalid data structure")
  if (is.matrix(x) | is.data.frame(x)) {
    if (ncol(x) != 1 & nrow(x) != 1)
      stop("invalid data structure")
  }
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid data type")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous data type")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  BootstrapFun.abun <- function(x, FunName, datatype, B) {
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2,
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- 1 - f1/n * A
    f0 <- max(round(f0.hat), 1)
    if (f0.hat == 0) {
      lambda <- 0
      if (sum(x > 0) == 1) {
        warning("The Bootstrap community has only one species. Estimation is not robust.")
      }
    }
    else {
      lambda <- (1 - Chat)/sum(x/n * (1 - x/n)^n)
    }
    pi <- x/n * (1 - lambda * (1 - x/n)^n)
    pi.star <- c(pi, rep((1 - Chat)/f0, f0))
    X <- rmultinom(B, n, pi.star)
    se <- sd(apply(X, 2, function(x) FunName(x, datatype)))
    return(se)
  }
  BootstrapFun.ince <- function(y, FunName, datatype, B) {
    t <- y[1]
    y <- y[-1]
    y <- y[y > 0]
    U <- sum(y)
    Q1 <- sum(y == 1)
    Q2 <- sum(y == 2)
    Q0.hat <- ifelse(Q2 == 0, (t - 1)/t * Q1 * (Q1 - 1)/2,
                     (t - 1)/t * Q1^2/2/Q2)
    A <- ifelse(Q1 > 0, t * Q0.hat/(t * Q0.hat + Q1), 1)
    Chat <- 1 - Q1/U * A
    Q0 <- max(round(Q0.hat), 1)
    if (Q0.hat == 0) {
      tau <- 0
      if (sum(y > 0) == 1) {
        warning("The Bootstrap community has only one species. Estimation is not robust.")
      }
    }
    else {
      tau <- U/t * (1 - Chat)/sum(y/t * (1 - y/t)^t)
    }
    pi <- y/t * (1 - tau * (1 - y/t)^t)
    pi.star <- c(pi, rep(U/t * (1 - Chat)/Q0, Q0))
    y1 <- rbind(t, matrix(rbinom(length(pi.star) * B, t,
                                 pi.star), ncol = B))
    tmp <- which(colSums(y1) == t)
    if (length(tmp) > 0)
      y1 <- y1[, -tmp]
    se <- sd(apply(y1, 2, function(y2) FunName(y2, datatype)))
    return(se)
  }
  if (datatype == "abundance") {
    BootstrapFun.abun(x = x, FunName, datatype, B)
  }
  else if (datatype == "incidence") {
    BootstrapFun.ince(y = x, FunName, datatype, B)
  }
}
DataInfo = function (x, datatype = "abundance") {
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  Fun.abun <- function(x) {
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x == k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x > 0)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2,
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- round(1 - f1/n * A, 4)
    c(n, Sobs, Chat, fk)
  }
  Fun.ince <- function(x) {
    nT <- x[1]
    x <- x[-1]
    U <- sum(x)
    Qk <- sapply(1:10, function(k) sum(x == k))
    Q1 <- Qk[1]
    Q2 <- Qk[2]
    Sobs <- sum(x > 0)
    Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2,
                     (nT - 1)/nT * Q1^2/2/Q2)
    A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1),
                1)
    Chat <- round(1 - Q1/U * A, 4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  if (datatype == "abundance") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.abun(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.abun))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.abun))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "n", "S.obs", "SC", paste("f",
                                                           1:10, sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "n", "S.obs", "SC", paste("f",
                                                           1:10, sep = ""))
    }
    as.data.frame(out)
  }
  else if (datatype == "incidence") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.ince(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.ince))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.ince))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "T", "U", "S.obs", "SC",
                         paste("Q", 1:10, sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "T", "U", "S.obs", "SC",
                         paste("Q", 1:10, sep = ""))
    }
    as.data.frame(out)
  }
}



ggiNEXT = function (x, type = 1, se = TRUE, facet.var = "none", color.var = "site",
                    grey = FALSE) {
  UseMethod("ggiNEXT", x)
}

ggiNEXT.iNEXT = function (x, type = 1, se = TRUE, facet.var = "none", color.var = "site",
                          grey = FALSE) {
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var,
                                               SPLIT) == -1)
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var,
                                               SPLIT) == -1)
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order")
    color.var <- "site"
  if (facet.var == "site")
    color.var <- "order"
  options(warn = -1)
  z <- fortify(x, type = type)
  options(warn = 0)
  if (ncol(z) == 7) {
    se <- FALSE
  }
  datatype <- unique(z$datatype)
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }
    else {
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$shape <- factor(z$order)
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep = "-")
  }
  zz = z
  z$method[z$method == "observed"] = "interpolated"
  z$lty <- z$lty <- factor(z$method, levels = unique(c("interpolated",
                                                       "extrapolated"), c("interpolation", "interpolation",
                                                                          "extrapolation")))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method == "observed"), ]
  g <- ggplot(z, aes_string(x = "x", y = "y", colour = "col")) +
    geom_point(aes_string(shape = "shape"), size = 5, data = data.sub)
  g <- g + geom_line(aes_string(linetype = "lty"), lwd = 1.5) +
    guides(linetype = guide_legend(title = "Method"), colour = guide_legend(title = "Guides"),
           fill = guide_legend(title = "Guides"), shape = guide_legend(title = "Guides")) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          text = element_text(size = 18), legend.key.width = unit(1.2,
                                                                  "cm"))
  if (type == 2L) {
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Sample coverage")
  }
  else if (type == 3L) {
    g <- g + labs(x = "Sample coverage", y = "Species diversity")
  }
  else {
    g <- g + labs(x = "Number of sampling units", y = "Species diversity")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Species diversity")
  }
  if (se)
    g <- g + geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr",
                                    fill = "factor(col)", colour = "NULL"), alpha = 0.2)
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1 & type != 2) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_wrap(~order, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              ncol = length(levels(factor(z$order))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_wrap(~site, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in%
        names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$site))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (grey) {
    g <- g + theme_bw(base_size = 18) + scale_fill_grey(start = 0,
                                                        end = 0.4) + scale_colour_grey(start = 0.2, end = 0.2) +
      guides(linetype = guide_legend(title = "Method"),
             colour = guide_legend(title = "Guides"), fill = guide_legend(title = "Guides"),
             shape = guide_legend(title = "Guides")) + theme(legend.position = "bottom",
                                                             legend.title = element_blank())
  }
  g <- g + theme(legend.box = "vertical")
  return(g)
}

ggiNEXT.default = function (x, ...) {
  stop("iNEXT doesn't know how to deal with data of class ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

DataInfo = function (x, datatype = "abundance") {
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "list") {
      x <- lapply(x, as.incfreq)
    }
    else {
      x <- as.incfreq(x)
    }
    datatype <- "incidence"
  }
  Fun.abun <- function(x) {
    n <- sum(x)
    fk <- sapply(1:10, function(k) sum(x == k))
    f1 <- fk[1]
    f2 <- fk[2]
    Sobs <- sum(x > 0)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2,
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- round(1 - f1/n * A, 4)
    c(n, Sobs, Chat, fk)
  }
  Fun.ince <- function(x) {
    nT <- x[1]
    x <- x[-1]
    U <- sum(x)
    Qk <- sapply(1:10, function(k) sum(x == k))
    Q1 <- Qk[1]
    Q2 <- Qk[2]
    Sobs <- sum(x > 0)
    Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2,
                     (nT - 1)/nT * Q1^2/2/Q2)
    A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1),
                1)
    Chat <- round(1 - Q1/U * A, 4)
    out <- c(nT, U, Sobs, Chat, Qk)
  }
  if (datatype == "abundance") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.abun(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.abun))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.abun))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "n", "S.obs", "SC", paste("f",
                                                           1:10, sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "n", "S.obs", "SC", paste("f",
                                                           1:10, sep = ""))
    }
    as.data.frame(out)
  }
  else if (datatype == "incidence") {
    if (class(x) == "numeric" | class(x) == "integer") {
      out <- matrix(Fun.ince(x), nrow = 1)
    }
    else if (class(x) == "list") {
      out <- do.call("rbind", lapply(x, Fun.ince))
    }
    else if (class(x) == "matrix" | class(x) == "data.frame") {
      out <- t(apply(as.matrix(x), 2, Fun.ince))
    }
    if (nrow(out) > 1) {
      out <- data.frame(site = rownames(out), out)
      colnames(out) <- c("site", "T", "U", "S.obs", "SC",
                         paste("Q", 1:10, sep = ""))
      rownames(out) <- NULL
    }
    else {
      out <- data.frame(site = "site.1", out)
      colnames(out) <- c("site", "T", "U", "S.obs", "SC",
                         paste("Q", 1:10, sep = ""))
    }
    as.data.frame(out)
  }
}

as.abucount = function (x) {
  if (class(x) == "data.frame" | class(x) == "matrix") {
    y <- rowSums(x)
    y
  }
  else if (class(x) == "numeric" | class(x) == "integer" |
           class(x) == "double") {
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    x
  }
  else {
    stop("invalid data type, it should be a data.frame or matrix.")
  }
}

as.incfreq = function (x) {
  if (class(x) == "data.frame" | class(x) == "matrix") {
    a <- sort(as.numeric(unique(c(unlist(x)))))
    if (!identical(a, c(0, 1))) {
      warning("Invalid data type, the element of species by sites presence-absence matrix should be 0 or 1. Set nonzero elements as 1.")
      x <- (x > 0)
    }
    nT <- ncol(x)
    y <- rowSums(x)
    y <- c(nT, y)
    y
  }
  else if (class(x) == "numeric" | class(x) == "integer" |
           class(x) == "double") {
    warnings("Ambiguous data type, the input object is a vector. Set total number of sampling units as 1.")
    c(1, x)
  }
  else {
    stop("Invalid data type, it should be a data.frame or matrix.")
  }
}

estimateD = function (x, datatype = "abundance", base = "size", level = NULL,
                      conf = 0.95) {
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if (is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  if (datatype == "incidence") {
    stop("datatype=\"incidence\" was no longer supported after v2.0.8, \n         please try datatype=\"incidence_freq\".")
  }
  if (datatype == "incidence_freq")
    datatype <- "incidence"
  if (datatype == "incidence_raw") {
    if (class(x) == "data.frame" | class(x) == "matrix")
      x <- as.incfreq(x)
    else if (class(x) == "list")
      x <- lapply(x, as.incfreq)
    datatype <- "incidence"
  }
  BASE <- c("size", "coverage")
  if (is.na(pmatch(base, BASE)))
    stop("invalid datatype")
  if (pmatch(base, BASE) == -1)
    stop("ambiguous datatype")
  base <- match.arg(base, BASE)
  if (base == "size") {
    tmp <- invSize(x, datatype, size = level, conf = conf)
  }
  else if (base == "coverage") {
    tmp <- invChat(x, datatype, C = level, conf = conf)
  }
  tmp <- tmp[!duplicated(tmp), ]
  nam <- names(x)
  if (is.null(nam)) {
    tmp
  }
  else if (ncol(tmp) == 6) {
    tmp <- cbind(site = nam, tmp)
  }
  else {
    tmp <- cbind(site = rep(nam, each = 3), tmp)
  }
  rownames(tmp) <- NULL
  tmp
}

fortify.iNEXT = function (model, data = model$iNextEst, type = 1, ...) {
  datatype <- ifelse(names(model$DataInfo)[2] == "n", "abundance",
                     "incidence")
  z <- data
  if (class(z) == "list") {
    z <- data.frame(do.call("rbind", z), site = rep(names(z),
                                                    sapply(z, nrow)))
    rownames(z) <- NULL
  }
  else {
    z$site <- ""
  }
  if (ncol(z) == 6) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }
  else if (ncol(z) > 6) {
    se <- TRUE
  }
  if (type == 1L) {
    z$x <- z[, 1]
    z$y <- z$qD
    if (se) {
      z$y.lwr <- z[, 5]
      z$y.upr <- z[, 6]
    }
  }
  else if (type == 2L) {
    if (length(unique(z$order)) > 1) {
      z <- subset(z, order == unique(z$order)[1])
    }
    z$x <- z[, 1]
    z$y <- z$SC
    if (se) {
      z$y.lwr <- z[, 8]
      z$y.upr <- z[, 9]
    }
  }
  else if (type == 3L) {
    z$x <- z$SC
    z$y <- z$qD
    if (se) {
      z$y.lwr <- z[, 5]
      z$y.upr <- z[, 6]
    }
  }
  z$datatype <- datatype
  z$plottype <- type
  if (se) {
    data <- z[, c("datatype", "plottype", "site", "method",
                  "order", "x", "y", "y.lwr", "y.upr")]
  }
  else {
    data <- z[, c("datatype", "plottype", "site", "method",
                  "order", "x", "y")]
  }
  data
}

print.iNEXT = function (x, ...) {
  site.n <- nrow(x$DataInfo)
  order.n <- ifelse(site.n > 1, paste(unique(x$iNextEst[[1]]$order),
                                      collapse = ", "), paste(unique(x$iNextEst$order), collapse = ", "))
  cat("Compare ", site.n, " assemblages with Hill number order q = ",
      order.n, ".\n", sep = "")
  cat("$class: iNEXT\n\n")
  cat("$DataInfo: basic data information\n")
  print(x$DataInfo)
  cat("\n")
  cat("$iNextEst: diversity estimates with rarefied and extrapolated samples.\n")
  if (class(x$iNextEst) == "data.frame") {
    y <- x$iNextEst
    m <- quantile(y[, 1], type = 1)
    res <- y[y[, 1] %in% m, ]
  }
  else {
    res <- lapply((x$iNextEst), function(y) {
      m <- quantile(y[, 1], type = 1)
      y[y[, 1] %in% m, ]
    })
  }
  print(res)
  cat("\n")
  cat("$AsyEst: asymptotic diversity estimates along with related statistics.\n")
  print(x$AsyEst)
  cat("\n")
  cat("NOTE: Only show five estimates, call iNEXT.object$iNextEst. to show complete output.\n")
  return(invisible())
}
