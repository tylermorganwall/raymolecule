#'@title Cross Product of two
#'
#'@param a Length-3 numeric vector.
#'@param b Length-3 numeric vector.
#'@keywords internal
cross = function(a, b) {
  i1 = c(2,3,1)
  i2 = c(3,1,2)
  return(a[i1]*b[i2] - a[i2]*b[i1])
}

#'@title Unit vector
#'
#'@param v Numeric vector.
#'@keywords internal
unit_vector = function(v) {
  return(v/sqrt(sum(v*v)))
}

#'@title Create Orthonormal Basis from w (z)
#'
#'@param v Numeric vector.
#'@keywords internal
onb_from_w = function(n) {
  a1 = unit_vector(n)
  a = c()
  if(a1[1] > 0.999999) {
    a = c(0,1,0)
  } else {
    a = c(0,0,1)
  }
  a2 = unit_vector(cross(a1,a))
  a3 = cross(a1, a2)
  return(matrix(c(a1,a2,a3),nrow=3,byrow=TRUE))
}
