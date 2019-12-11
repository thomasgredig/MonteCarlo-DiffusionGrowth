# Height-Height Correlation Function (published in 2013)
#
# Thomas Gredig et al, Journal of Physics: 
# Conference Series 417 (2013) 012069
#################################################

hhcf <- function(m, NUM.ITER, bins) {
  dimx = ncol(m)
  dimy = nrow(m)
  px1 = round(runif(NUM.ITER, min=1, max=dimx))
  py1 = round(runif(NUM.ITER, min=1, max=dimy))
  px2 = round(runif(NUM.ITER, min=1, max=dimx))
  py2 = round(runif(NUM.ITER, min=1, max=dimy))
  
  r = sqrt( (px2-px1)^2 + (py2-py1)^2 )
  g=rep(0,NUM.ITER)
  for(i in 1:NUM.ITER) {
    g[i] = (m[px1[i],py1[i]] - m[px2[i],py2[i]])^2
  }
  data.frame(r,g)
}