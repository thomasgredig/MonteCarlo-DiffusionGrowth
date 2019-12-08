#########################################
# growth of molecules on surface by diffusion
#
# (c) 2019 Thomas Gredig
#########################################

#  Model Parameters
##########################
library(ggplot2)

# Parameters
############
N = 40           # array size
J = 1           # interaction strength
conv.eq = 500   # convergence to equilibrium
conv = 500      # measurements
reInit = FALSE  # re-initialize for new temperature
path.FIGS = 'images'
path.DATA = 'data'


# Monte Carlo Step
mcstep <- function() {
  
}


# Array Initialization
######################
mol = matrix(data=sign(runif(N*N)-0.5), nrow=N)

# Sample Output
###############
rasterGraph(spin)
ggsave(file.path(path.FIGS,paste0('Ising1D-',N,'x',N,'-Random.png')), width=4,height=4,dpi=220)
computeIsing1DRand(conv.eq*N*N, J, 1/2.2)
totalEnergy(N,J)
rasterGraph(spin)
ggsave(file.path(path.FIGS,paste0('Ising1D-',N,'x',N,'-Domains.png')), width=4,height=4,dpi=220)

# Computation Intesive Part: M vs T
##################################
d.runTime = data.frame(N,J,conv.eq,conv, reInit, date=Sys.Date(), 
                       start.time=as.numeric(Sys.time()),
                       end.time=0,diff.s=0)
bSeq = 1/TSeq
result = data.frame()
for(b in bSeq) {
  print(paste("Temp: ",1/b))
  if (reInit) { spin = matrix(data=sign(runif(N*N)-0.5), nrow=N) }
  computeIsing1DRand(conv.eq*N*N, J, b)
  Mavg = 0
  M2avg = 0
  for(i in 1:conv) {
    computeIsing1DRand(N*N, J, b)
    Ms = sum(spin)
    Mavg = Mavg + Ms
    M2avg = M2avg + Ms*Ms
  }
  result = rbind(result, data.frame(b, conv, Mavg, M2avg,
                                    T.J = 1/b/J,
                                    chi=(M2avg/(conv*N*N) - Mavg*Mavg/(conv*conv*N*N))*b))
}

# save timing
d.runTime$end.time = as.numeric(Sys.time())
d.runTime$diff.s = d.runTime$end.time-d.runTime$start.time
d.runTime$type = '1D'
d.runTime$step = TSeq[2]-TSeq[1]
d.runTime$TempStart = TSeq[1]
d.runTimeAll = rbind(d.runTimeAll,d.runTime)
write.csv(d.runTimeAll,file=file.runTime,row.names = FALSE)

# Graphing of Data
##################
ggplot(result, aes(T.J, abs(Mavg)/(conv*N*N))) +
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab('|M|') +
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising1D-',N,'x',N,'-c',conv,'.png')), width=4, height=3, dpi=220)
write.csv(result,file.path(path.DATA,paste0('Ising1D-',N,'x',N,'-c',conv,'.csv')), row.names=FALSE)

ggplot(result, aes(T.J, chi)) +
  geom_smooth(span=0.2)+
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab(expression(paste(chi))) +
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising1D-',N,'x',N,'-c',conv,'-Chi.png')), width=4, height=3, dpi=220)

