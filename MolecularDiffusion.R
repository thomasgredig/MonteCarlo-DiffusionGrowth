#########################################
# growth of molecules on surface by diffusion
#
# (c) 2019 Thomas Gredig
#########################################

#  Model Parameters
##########################
library(ggplot2)
library(cowplot)

# Parameters
############
N = 500           # array size
diffusionSteps = 1000
conv.eq = 500   # convergence to equilibrium
conv = 500      # measurements
reInit = FALSE  # re-initialize for new temperature
path.FIGS = 'images'
path.DATA = 'data'


# Monte Carlo Step
molDiffusion <- function(x,y) {
  steps=0
  while(steps<diffusionSteps) {
    nb=mol[(x %% N)+1,y] + mol[((x-2) %% N)+1,y] + 
      mol[x,(y %% N)+1] + mol[x,((y-2) %% N)+1] 
    if(nb>0) break
    # otherwise do a diffusion step
    if (runif(1)>0.5) { 
      x = ((x+round(runif(1,0,1))*2) %% N) + 1
    } 
    if (runif(1)>0.5) { #else {
      y = ((y+round(runif(1,0,1))*2 ) %% N) + 1
    }
    steps=steps+1
  }
  mol[x,y] <<- 1
}

# adding one molecule
molAdd <- function() {
  if (sum(mol)>N*N) break
  while(1) {
    x=round(runif(1,min=1,max=N))
    y=round(runif(1,min=1,max=N))
    if(mol[x,y]==0) break
  }
  molDiffusion(x,y)
}

rasterGraph <- function(N, spinMatrix) {
  df = expand.grid(x = 1:N, y = 1:N)
  df$spin = as.vector(spinMatrix)
  ggplot(df, aes(x,y, fill=spin)) + geom_raster() + 
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colors=c("red","blue"))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw() + theme(legend.position='none')
}



# Array Initialization
######################
mol = matrix(data=0, nrow=N, ncol=N)
for(i in 1:(N/2)) {
  for(j in 1:N) {
    molAdd()
  }
  # Sample Output
  ###############
  #print(rasterGraph(N,mol))
  #sum(mol)
  #Sys.sleep(0.1)
}


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


