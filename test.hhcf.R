# Testing the HHCF
library(ggplot2)
library(cowplot)
source('hhcf.func.R')
source('func.R')
path.FIGS = 'images'
# Parameters
############
N = 400           # array size
diffusionSteps = 1000

# Array Initialization
######################
mol = matrix(data=0, nrow=N, ncol=N)
for(i in 1:(N/2)) {
  for(j in 1:N) {
    molAdd()
  }
}

NUM.ITER = 1e5
#r = data.frame()
for(NUM.ITER in c(5:100*1e5)) {
  q = hhcf(mol, NUM.ITER)
  bs = seq(0, 120, length.out = 500)
  q$bin = cut(q$r, breaks=bs)
  bg = as.vector(by(q$g, q$bin, sum, simplify = TRUE))
  bg.norm = as.vector(by(q$g, q$bin, length, simplify = TRUE))
  bg1 = bg/bg.norm
  df = data.frame(
    r = bs,
    g = c(0,bg1)
  )
  df1 = subset(df, r<max(bs)*0.4)
  
  nls(data=df1,
      g ~ A*(1-exp(-(r/xi))),
      start=list(A=0.45, xi=5))->fit
  df1$g.fit = predict(fit, list(r=df1$r))
  xi = signif(summary(fit)$coeff[2],3)
  xi.sd = signif(summary(fit)$coeff[4],3)
  sigma= signif(sqrt(summary(fit)$coeff[1]),3)
  r = rbind(r, data.frame(N, diffusionSteps, num = NUM.ITER, sigma, xi,xi.sd))
  print(paste0("iter=", NUM.ITER))
}
plot(r$num, r$xi)

r$diffusionSteps = factor(r$diffusionSteps)
write.csv(r, file.path('data','test.hhcf.results.csv'), row.names=FALSE)

ggplot(r, aes(num/1e6, xi, color=diffusionSteps)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=xi-xi.sd, ymax=xi+xi.sd)) + 
  scale_x_continuous(breaks=0:5*2) + 
  ylab(expression(paste(xi))) +
  theme_bw() +
  xlab(expression(paste('number of iterations (10'^6,')')))
ggsave(file.path(path.FIGS,paste0('test.hhcf-',N,'x',N,'-HHCF',diffusionSteps,'.png')), 
       width=6, height=4, dpi=300)

ggplot(r, aes(num/1e6, xi.sd/xi*100, color=diffusionSteps)) + 
  geom_point() + theme_bw() +
  scale_y_continuous(limits=c(0,13))
