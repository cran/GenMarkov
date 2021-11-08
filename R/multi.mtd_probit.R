
#MMC estimation with Nicolau (2014) specification through maxLik()
multi.mtd_probit <- function(y, initial, nummethod = "bfgs") {
  results <- list(NA)


  #Calculate Frequency Transition Matrices (Ching (2002))
  ProbMatrixF <- function(s) {
    n = nrow(s)
    m1 = max(s)
    m0 = min(s)
    f = array(1, dim = c(m1, m1, (ncol(s) * ncol(s))))
    d = 0

    for (i in 1:ncol(s)) {
      for (j in 1:ncol(s)) {
        d = d + 1
        for (k1 in m0:m1) {
          for (k2 in m0:m1) {
            c = 0
            for (t in 2:n) {
              if (s[t - 1, j] == k1 && s[t, i] == k2) {
                c = c + 1
              }
            }
            f[k2, k1, d] = c
          }
        }
      }
    }
    return(f)
  }

  #Calculate Probability Transition Matrices (Ching (2002))
  ProbMatrixQ <- function(f, s) {
    n = nrow(s)
    m1 = max(s)
    m0 = min(s)
    q = array(0, dim = c(m1, m1, (ncol(s) * ncol(s))))
    d = 0

    for (i in 1:ncol(s)) {
      for (j in 1:ncol(s)) {
        d = d + 1
        for (k1 in m0:m1) {
          for (k2 in m0:m1) {
            if (sum(f[, k1, d]) > 0) {
              q[k2, k1, d] = f[k2, k1, d] / sum(f[, k1, d])
            }
          }
        }
      }
    }

    return(q)
  }

  #Array of number of sequences: n_i0_il
  NumberOfSequences <- function(f, s) {
    m1 <- max(s)

    n = m1 * m1

    n_i0_il <- array(0, dim = c(n, (ncol(s) * ncol(s))))

    for (g in 1:(ncol(s) * ncol(s))) {
      n_i0_il[, g] = matrixcalc::vec(f[, , g]) #Vectorization of frequency transition matrices
    }
    return(n_i0_il)
  }

  #Array of the probabilities values: q_i0_il
  ArrayQ <- function(q, s) {
    m1 <- max(s)
    n = m1 * m1
    q_i0_il <- array(0, dim = c(n, (ncol(s) * ncol(s))))

    for (g in 1:(ncol(s) * ncol(s))) {
      q_i0_il[, g] = matrixcalc::vec(q[, , g]) #Vectorization of probability transition matrices
    }
    return(q_i0_il)
  }


  nr.eq <- ncol(y)
  m1 <- max(y)
  ini <- initial

  f <- ProbMatrixF(y)
  q <- ProbMatrixQ(f, y)

  ni <- NumberOfSequences(f, y)
  qi <- ArrayQ(q, y)

  j = 0
  k = 0
  for (i in 1:nr.eq) {
    ll <- 0

    neq <- ni[, (i + j):(nr.eq * i)]
    qeq <- qi[, (i + j):(nr.eq * i)]

    #Define Log-likelihood as in Nicolau (2014)
    LogLikelihood_p <-
      function(eta,
               n_i0_il = ni,
               q_i0_il = qi,
               n = neq,
               qn = qeq) {
        nr.eq <- length(eta) - 1
        p <- array(1, dim = c(1, nr.eq))

        k = 0
        for (i in 1:nr.eq) {
          p[i] <-
            stats::pnorm(sum(cbind(rep(
              1, length(q_i0_il[, i])
            ), q_i0_il[, (i + k):(nr.eq * i)]) %*% eta))
          k = k + nr.eq - 1
        }

        ll <- 0
        for (j in 1:ncol(n)) {
          for (i in 1:length(n[, j])) {
            ll <-
              ll + n[i, j] * log(stats::pnorm(append(qn[i, ], 1, after = 0) %*% eta) / sum(p))
          }
        }

        ll
      }


    otim <- maxLik::maxLik(LogLikelihood_p, start = ini, method = nummethod)


    res <- summary(otim)

    results[[i + k]] <- res$estimate
    results[[(i + k + 1)]] <- res$loglik
    names(results)[[i + k]] <- paste("Equation", i, sep = " ")
    names(results)[[(i + k + 1)]] <- paste("LogLik", i, sep = " ")

    j = j + nr.eq - 1
    k = k + 1
  }

  return(results)

}
