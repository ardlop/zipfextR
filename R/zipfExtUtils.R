# Hurwitz Zeta function ------------------------------

.zeta_x<-function(alpha, x) {
  if(!is.numeric(alpha) || !is.numeric(x)){
    stop("Wrong input values!")
  }

  if(x < 1){
    stop("Error!! The 'x' have to be greater or iqual than one.")
  }

  aux <- 0
  if(x == 1) {
    VGAM::zeta(alpha)
  }
  else {
    VGAM::zeta(alpha) - sum((1:(x-1))^(-alpha))
  }
}
