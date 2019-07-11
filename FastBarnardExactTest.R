########################################################################
# fast algorithm for Barnard Exact Test code
# This is efficientised from the source code in
# https://www.r-bloggers.com/barnard%E2%80%99s-exact-test-%E2%80%93-a-powerful-alternative-for-fisher%E2%80%99s-exact-test-implemented-in-r/
# Please feel free to reuse the code with mentioning the source (Hamed Haselimashhadi, https://github.com/Hamedhm/Scatter-R-codes)
########################################################################

fastBarnardextest  = function(x,
                              tail = 2 ,
                              prob = seq(10 ^ -10, 1 - 10 ^ -10, length.out = 101),
                              plot = TRUE) {
  if (is.null(x)  ||
      !(is.table(x) ||
        is.matrix(x)) ||
      any(dim(x) != 2))
    stop('The input must be a 2 by 2 table/matrix')

  fprob = function(i, j, c1, c2) {
    n  = c1 + c2
    pa = i / c1
    pb = j / c2
    px = (i + j) / n
    if (px == 0 || pa == pb) {
      return(0)
    } else
      return((pa - pb) / sqrt(px * (1 - px) * ((1 / c1) + (1 / c2))))
  }
  c  = colSums(x)
  r  = rowSums(x)
  c1 = c[1]
  c2 = c[2]
  n = sum(c)
  pao = x[1, 1] / c1
  pbo = x[1, 2] / c2
  pxo = r[1] / n
  TXO = abs(pao - pbo) / sqrt(pxo * (1 - pxo) * (1 / c1 + 1 / c2))
  cbn = matrix(c(rep(0:c1, each = c2 + 1), rep(0:c2, c1 + 1)), nrow = 2, byrow = TRUE)

  n1     = lfactorial(c1)
  n2     = lfactorial(c2)
  lprob  = log(prob)
  clprob = log(1 - prob)
  Fact   = lfactorial(0:max(c1, c2, na.rm = TRUE)[1])
  ###################################################
  T = sapply(1:ncol(cbn), function(col) {
    i   = cbn[1, col]
    j   = cbn[2, col]
    s = n1 + n2 +
      (i + j) * lprob +
      (n - (i + j)) * clprob - sum(Fact[c(i + 1,
                                          j + 1,
                                          c1 - i + 1,
                                          c2 - j + 1)])
    t = fprob(i = i,
              j = j,
              c1 = c1,
              c2 = c2)
    return(c(t = t, s = exp(s)))
  })
  r = t(cbind(P = apply(T[-1,], 1, function(x) {
    sum(x[T[1,] >= TXO])
  }), prob = prob))
  ###################################################
  Nuisance.parameter = unlist(r[2, ][which.max(r[1, ])])
  p.value            = unlist(r[1, ][which.max(r[1, ])])
  if (plot) {
    plot(
      r[2, ],
      r[1, ],
      type = "l",
      main = "Barnard's exact P-value",
      xlab = "Nuisance parameter",
      ylab = "P.value"
    )
    abline(v = Nuisance.parameter, col = 2, lwd = 2)
  }
  return(
    list(
      p.value            = unname(min(tail * p.value, 1)),
      Nuisance.parameter = unname(Nuisance.parameter),
      Wald.Statistic     = unname(TXO)               ,
      tail               = tail                      ,
      seq                = r
    )
  )
}

##########################
### Example
x = matrix(c(
  a = 85,
  b = 50,
  c = 30,
  d = 25
), 2, 2, byrow = 1)


b = fastBarnardextest (x, plot = TRUE)
f = fisher.test(x)
###
b$p.value
f
