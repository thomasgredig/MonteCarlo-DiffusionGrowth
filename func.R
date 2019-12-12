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

# Monte Carlo Step
molDiffusionLeftRight <- function(x,y) {
  steps=0
  while(steps<diffusionSteps) {
    nb=mol[(x %% N)+1,y] + mol[((x-2) %% N)+1,y] + 
      mol[x,(y %% N)+1] + mol[x,((y-2) %% N)+1] 
    if(nb>0) break
    # otherwise do a diffusion step
    if (runif(1)>0.5) { 
      x = ((x+round(runif(1,0,1))*2) %% N) + 1
    } else {
      y = ((y+round(runif(1,0,1))*2 ) %% N) + 1
    }
    steps=steps+1
  }
  mol[x,y] <<- 1
}

# adding one molecule
molAdd <- function(typeLeftRight = TRUE) {
  if (sum(mol)>N*N) break
  while(1) {
    x=round(runif(1,min=1,max=N))
    y=round(runif(1,min=1,max=N))
    if(mol[x,y]==0) break
  }
  if (typeLeftRight) { molDiffusion(x,y) } else { molDiffusionLeftRight(x,y) }
}


# rastering graph
rasterGraph <- function(N, spinMatrix) {
  df = expand.grid(x = 1:N, y = 1:N)
  df$spin = as.vector(spinMatrix)
  ggplot(df, aes(x,y, fill=spin)) + geom_raster() + 
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colors=c("red","blue"))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw() + theme(legend.position='none')
}
