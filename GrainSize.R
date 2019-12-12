# Grain Size Dependence on Diffusion Length
library(ggplot2)
library(cowplot)
source('hhcf.func.R')
source('func.R')
path.FIGS = 'images'
path.RESULTS = 'data'

# Parameters
############
N = 100           # array size
diffusionSteps.list = c(10,seq(from=50, to=1000, length.out=20))

r = data.frame()
diffusionSteps = 50
for(diffusionSteps in diffusionSteps.list) {
  # Array Initialization
  ######################
  mol = matrix(data=0, nrow=N, ncol=N)
  for(i in 1:(N/2)) {
    for(j in 1:N) {
      molAdd(typeLeftRight=FALSE)
    }
  }

  NUM.ITER = 6e6
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
      start=list(A=0.667, xi=3))->fit
  df1$g.fit = predict(fit, list(r=df1$r))
  xi = signif(summary(fit)$coeff[2],3)
  xi.sd = signif(summary(fit)$coeff[4],3)
  sigma= signif(sqrt(summary(fit)$coeff[1]),3)
  r = rbind(r, data.frame(N, diffusionSteps, NUM.ITER, sigma, xi,xi.sd))
  print(paste0("diffusionSteps=", diffusionSteps))
}
plot(r$diffusionSteps, r$xi)

write.csv(r, file.path(path.RESULTS,'GrainSize.results.lf.csv'), row.names=FALSE)

ggplot(r, aes(sqrt(diffusionSteps),xi)) +
  geom_point(size=3) +
  geom_point(size=2, col='red') +
  scale_x_continuous(limits=c(0,35)) + 
  scale_y_continuous(limits=c(0,7)) + 
  ylab(expression(paste(xi))) +
  theme_bw() +
  xlab(expression(paste('square root of diffusion steps (D'^0.5,')')))
ggsave(file.path(path.FIGS,paste0('GrainSize-',N,'x',N,'-HHCF',diffusionSteps,'.png')),
       width=4, height=3, dpi=300)

