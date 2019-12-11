#########################################
# growth of molecules on surface by diffusion
#
# (c) 2019 Thomas Gredig
#########################################

#  Model Parameters
##########################
library(ggplot2)
library(cowplot)
source('hhcf.func.R')
source('func.R')

# Parameters
############
N = 400           # array size
diffusionSteps = 2000
reInit = FALSE  # re-initialize for new temperature
path.FIGS = 'images'
path.DATA = 'data'


# Array Initialization
######################
mol = matrix(data=0, nrow=N, ncol=N)
for(i in 1:(N/2)) {
  for(j in 1:N) {
    molAdd()
  }
}

q = hhcf(mol, 1e7)
bs = seq(0, 120, length.out = 500)
q$bin = cut(q$r, breaks=bs)
bg = as.vector(by(q$g, q$bin, sum, simplify = TRUE))
bg.norm = as.vector(by(q$g, q$bin, length, simplify = TRUE))
bg1 = bg/bg.norm
#plot(bs,c(0,bg1), xlim=c(0,max(bs)/2), main='Correlation Length')
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
sigma= signif(sqrt(summary(fit)$coeff[1]),3)

ggplot(df1, aes(r, g)) + 
  geom_point(col='black',size=2) + 
  geom_point(col='red') + 
#  scale_x_log10() + 
  geom_path(aes(r,g.fit), col='blue') + 
  ylab('g(r)') +
  ggtitle(paste0('FIT for ',N,'x',N,' (diff=',diffusionSteps,')'),
    subtitle = paste0('s=',sigma,' xi=',xi)) + 
  theme_bw()

ggsave(file.path(path.FIGS,paste0('MolDiff-Side-',N,'x',N,'-D',diffusionSteps,'-CorrLen.png')), 
       width=4, height=3, dpi=300)


theme_georgia <- function(...) {
  theme_gray(base_family = "Georgia", ...) + 
    theme(plot.title = element_text(face = "plain"))
}

g1 = rasterGraph(N,mol)
title_theme <- ggdraw() +
  draw_label(paste("Molecule Diffusion Simulation : \nN:",N,'x',N,"  Diffusion: ",diffusionSteps), 
             fontfamily = theme_georgia()$text$family, 
             fontface = theme_georgia()$plot.title$face, x = 0.05, hjust = 0)
  
plot_grid(title_theme,g1, ncol = 1, rel_heights = c(0.15, 1))
ggsave(file.path(path.FIGS,paste0('MolDiff-Side-',N,'x',N,'-D',diffusionSteps,'.png')), 
       width=4, height=4, dpi=300)


