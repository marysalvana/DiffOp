calculateDistanceC <- function (y, x) {
  n <- length(y)
  dist <- rep (0, n*n)
  res <- .C("calculatedistance", y = as.double(y), x = as.double(x), dist = as.double(dist), n = as.integer(n))
  dist <- matrix(res$dist, n, n)
  return (dist)
}

calculateDistanceEarth <- function (y, x) {
  n <- length(y)
  dist <- rep (0, n*n)
  res <- .C("DistanceEarth", y = as.double(y), x = as.double(x),
            dist = as.double(dist), n = as.integer(n))
  dist <- matrix(res$dist, n, n)
  return (dist)
}

CalculateCovBiDifferential <- function (x, y, z, a1, b1, c1, d1, a2, b2, c2, d2,
                                        scale_h, scale_v, nu, sigsq1, sigsq2, beta,
                                        splines_degree, inner_knots1, inner_knots2) {

  n <- length(y)

  dist <- rep(0, (2*n)*(2*n))

  if(splines_degree > 0){
    basis1 <- bsplineBasis(z, splines_degree, inner_knots1)
    nb1 <- ncol(basis1)
    basis2 <- bsplineBasis(z, splines_degree, inner_knots2)
    nb2 <- ncol(basis2)

    c1 <- c(basis1 %*% matrix(c1, ncol = 1))
    c2 <- c(basis2 %*% matrix(c2, ncol = 1))
  }else{
    c1 <- rep(c1, n)
    c2 <- rep(c2, n)
  }

  res <- .C("CovBiDifferential", x = as.double(x), y = as.double(y), z = as.double(z),
            a1 = as.double(a1), b1 = as.double(b1), c1 = as.double(c1), d1 = as.double(d1),
            a2 = as.double(a2), b2 = as.double(b2), c2 = as.double(c2), d2 = as.double(d2),
            scale_h = as.double(scale_h), scale_v = as.double(scale_v), nu = as.double(nu),
            sigsq1 = as.double(sigsq1), sigsq2 = as.double(sigsq2), beta = as.double(beta),
            dist = as.double(dist), n = as.integer(n))

  dist <- matrix(res$dist, 2*n, 2*n)
  return (dist)
}
