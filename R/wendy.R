# setwd(directory)
# X <- read.csv("testdata.csv", header=TRUE,skipNul=TRUE) # Read from file

wendy <- function(file, ...) {
  library(deSolve) # has to be installed before first use
  library(ggplot2) # has to be installed before first use
  if(!grepl(".csv$",file)) { # check if the file is csv file
    stop("Can only upload .csv file!")
  }
  X <- read.csv(file, ...) # Read from file
  H<-X$N_FORAGES[1] # H = number of forages = number of habitats
  forages <-as.character(X$FORAGE_NAME[1:H]) # Names of energy (food) types
  Tarea<-X$TOTAL_AREA_BY_FORAGE[1:H] # Total area for each habitat or forage
  AvtimeColumn <- grep("AVAIL_TIME",colnames(X)) # time column for availability
  Availtime<- X[,AvtimeColumn][!is.na(as.numeric(X[,AvtimeColumn]))]
  Habav <- X[1:length(Availtime),(AvtimeColumn+1):(AvtimeColumn+H)] # Availability
  I<-X$N_FOODTYPES[1] # I = number of energy types
  foods <-as.character(X$FOOD_TYPE[1:I]) # Names of energy (food) types
  DensityColumn <- grep("FOOD_TYPE",colnames(X))+1 # Energy density must be next to foods column
  Den<-X[1:I,DensityColumn:(DensityColumn+H-1)] # read energy density by food type for each habitat
  K<-X$N_GUILDS[1] # K = number of water fowl guilds
  guilds <-as.character(X$GUILD_NAME[1:K]) # Names of waterfowl guilds
  PrefColumn<-grep("GUILD_NAME",colnames(X))+1 # Preference info must be next to guilds column
  Pref<-X[1:I,PrefColumn:(PrefColumn+K-1)] # Preference matrix
  PoptimeColumn <- grep("POP_TIME",colnames(X)) # time column for population
  Poptime <- X[,PoptimeColumn][!is.na(as.numeric(X[,PoptimeColumn]))]
  # as.numeric(gsub(",", "", trades1$Price))
  Pop<-X[1:length(Poptime),(PoptimeColumn+1):(PoptimeColumn+K)] # guild popolution
  IntaketimeColumn <- grep("INTAKE_TIME",colnames(X)) # time column for daily intake by guild
  Intaketime <- X[,IntaketimeColumn][!is.na(as.numeric(X[,IntaketimeColumn]))]
  Intake <- X[1:length(Intaketime),(IntaketimeColumn+1):(IntaketimeColumn+K)]
  Ek <-as.numeric(Intake[1,1:K])
  Pars <- c(d<- 0.001,b<- 0.001) # Set decomposition rates
  # correction for input error when availability exceeds maximum area
  for (h in 1:H) {
    maxavail=max(Habav[,h]) # max availability
    if (maxavail>Tarea[h]) {
      print("Availability exceeds maximum area. Habitat area is reset")
      Tarea[h]<-maxavail}
  }
  # plot the availability
  matplot(Availtime,Habav,type = c("b"),pch=1,col = 1:H,xlab='Days',ylab='Available Acres')
  legend("topright", legend = forages[1:H], col=1:H,pch=1) # optional legend
  # calculate the curve fit to the population and plot
  mdays=max(Poptime)
  dayz=1:mdays
  popplot<-matrix(nrow=length(dayz),ncol=K)
  for(daycount in 1:mdays)
    for (k in 1:K) {
      popul<-splinefun(Poptime,Pop[,k],method = "natural") # spline fit to population
      popplot[daycount,k]=max(c(popul(daycount,deriv=0),0)) # popul must be non-negative
    }
  # optional non-stacked population plot
  # matplot(dayz,popplot,type = c("b"),pch=19,col = 1:K,xlab='Days',ylab='Population')
  # legend("topright", legend = 1:K, col=1:K,pch=19) # non-stacked plot
  population=matrix(t(popplot),nrow=(K*mdays),ncol=1)
  da <- data.frame(day=rep(1:mdays,each=K),guild=rep(guilds),val=population)
  p<- ggplot(da, aes(day,population))
  p+geom_area(aes(colour=guild,fill=guild),position='stack')
  # end of plottting population
  # Initialization of matrices used in ODE definition
  c<-array(0,dim=c(H,I,K)) # create array of consumption matrix
  dE<-matrix(0,nrow=I,ncol=H) # create matrix of time derivative of E (energy)
  dR<-matrix(0,nrow=I,ncol=H) # create matrix of time derivative of R (reserve)
  E<-matrix(0,nrow=I,ncol=H)
  R<-matrix(0,nrow=I,ncol=H)
  n_birds<- matrix(0,nrow=K,ncol=1)
  ### TOP of ODE DEFINITION ###
  Ened <- function(Time, y, Pars){
    with(as.list(c(y,Pars)),{
      #Redistribution of initial data
      E<-matrix(y[1:(H*I)],nrow=I,ncol=H)
      R<-matrix(y[(H*I+1):(2*H*I)],nrow=I,ncol=H)
      #Determination of energy consumption coefficient
      pEkh=t(Pref) %*% E #matrix multiplication Transpose(p[i,k])*E[i,h]
      for(k in 1:K){
        n_birds[k]= approx(dayz,popplot[,k],xout=Time,method="linear")$y
        pEk=sum(pEkh[k,])+0.1
        intake<-splinefun(Intaketime,Intake[,k])
        Ek[k]<-intake(Time,deriv=0)
        #  if (pEk<0.1) { print("Energy exhaustion for the guild #, at time, ")
        #    print(c(k,Time,pEk))
        #    stop
        #  }
        for(h in 1:H){
          for(i in 1:I){
            c[h,i,k]<- Pref[i,k]*E[i,h]/pEk*Ek[k]
          }
        }
      }
      #Determine energy change between available and reserve energy
      for(h in 1:H){
        int<-findInterval(Time,Availtime) # find interval for linear interpolation
        dA_dt<- (Habav[int+1,h]-Habav[int,h])/(Availtime[int+1]-Availtime[int])
        A <- Habav[int,h]+dA_dt*(Time-Availtime[int])
        if(A<0.1) dA_dt=max(c(dA_dt,0)) # Area cannot decrease if it is zero
        if(Tarea[h]-A<0.1) dA_dt=min(c(dA_dt,0)) # Area cannot increase if near max
        for(i in 1:I){
          if(dA_dt>0.01){
            rhi<- R[i,h]*dA_dt/(Tarea[h]-A+.1)
          }else if(dA_dt> -0.01){
            rhi <- 0
          }else rhi<- E[i,h]/(A+0.1)*dA_dt # dA is negative, rhi is negative
          cPk<-sum(c[h,i,]*n_birds)
          dE[i,h]= -d * E[i,h]  + rhi  -cPk
          dR[i,h]= -b * R[i,h]- rhi
          if(i==I & foods[I]=="BENTHIC") #Foodtype I = invertebrate;
          {
            gammah <- Den[I,h] # carrying capacity (density) stored in Den[I,h]
            if (gammah<0.01) # negligible invertebrates
            {dE[I,h]=0
            dR[I,h]=0}
            else {
              deltah <- gammah # carrying capacity in Reservoir set equal to that of Available
              if (A>0.1) dE[I,h]=0.1*(1-E[I,h]/(A*gammah+0.1))*E[I,h]+rhi-cPk
              else dE[I,h]=rhi-cPk #no breeding if area =0
              if ((Tarea[h]-A)>0.1) dR[I,h]=0.1*(1-R[I,h]/((Tarea[h]-A)*deltah+0.1))*R[I,h]-rhi
              else dR[I,h]=-rhi #no breeding if area =0
            }
          }
        }
      }
      list(c(as.vector(dE),as.vector(dR)))
    })
  }
  # End of ODE definition
  #Calculation of initial energy available and in reserve
  for(h in 1:H){
    E[,h]<-Den[,h]*Habav[1,h] # row 1 of Habav is the initial available area
    R[,h]<-Den[,h]*(Tarea[h]-Habav[1,h])
  }
  y = c(as.vector(E),as.vector(R)) #Vectorization of matrices before passing to ODEs
  totalsteps= min(c(max(Availtime),max(Poptime),max(Intaketime)))-2
  # simu time determined by available data
  tstep <- 1 # time step is set to one day
  # Initialization of matrices for output
  days=matrix(0,nrow=totalsteps+1,ncol=1)
  EnergyAvailable=matrix(0,nrow=totalsteps+1,ncol=H*I)
  E_food_type=matrix(0,nrow=totalsteps+1,ncol=I)
  E_habitat=matrix(0,nrow=totalsteps+1,ncol=H)
  EnergyReserve=matrix(0,nrow=totalsteps+1,ncol=H*I)
  TotalDemand <- matrix(0,nrow=totalsteps+1,ncol=1)
  dailyDemand <- matrix(0,nrow=totalsteps+1,ncol=1)
  energydeficit <- matrix(0,nrow=totalsteps+1,ncol=1)
  EnergyAvailable[1,]=y[1:(H*I)] # record the first output: given by initial conditions
  E_food_type[1,]=rowSums(E)
  E_habitat[1,]=colSums(E)
  EnergyReserve[1,]=y[(H*I+1):(2*H*I)]
  for (istep in 1:totalsteps) {
    day <- seq(istep*tstep,(istep+1)*tstep)
    # check if energy is sufficient: compare available with dailyneed
    allavailable<-sum(y[1:(H*I)])
    for(k in 1:K){
      n_birds[k]= approx(dayz,popplot[,k],xout=(istep+1)*tstep,method="linear")$y
      Ek[k]<- approx(Intaketime,Intake[,k],xout=(istep+1)*tstep,method="linear")$y
    }
    dailyneed<-sum(Ek*n_birds)
    if (dailyneed>allavailable/1.5){
      print('Energy is insufficient. It has been increased to sustain the population')
      print(c(dailyneed,allavailable))
      rel_importance=rowSums(Pref)/sum(Pref) # relative importance determined by preference matrix
      y[1:(H*I)]<-dailyneed*rel_importance/H +y[1:(H*I)] # distribute added energy accordingly
      energydeficit[istep] <- dailyneed # track the added energy as deficit
    }
    # call to the ODE solver
    out <- ode(y, func=Ened, parms=Pars, times = day, method="rk4")
    len<-length(out[,1])
    days[istep+1]<-out[len,1] # save current time
    y=out[len,2:(2*H*I+1)] # update initial condition
    E<-matrix(y[1:(H*I)],nrow=I,ncol=H)
    R<-matrix(y[(H*I+1):(2*H*I)],nrow=I,ncol=H)
    EnergyAvailable[istep+1,]=y[1:(H*I)] # Remaining energies (I x H)
    E_food_type[istep+1,]=rowSums(E) # sum energy by food type
    E_habitat[istep+1,]=colSums(E) # sum energy by habitat
    EnergyReserve[istep+1,]=y[(H*I+1):(2*H*I)] #
    dailyDemand[istep+1] <-dailyneed
    TotalDemand[istep+1] <-TotalDemand[istep]+dailyneed
    if(any(EnergyAvailable<0)) print("At least one energy type exhausted")
    if(sum(y[1:(H*I)])<0) stop("All energy exhausted")
  }
  TotalAvailEnergy <- rowSums(EnergyAvailable)
  TotalResEnergy <- rowSums(EnergyReserve)
  TotalEnergy <- TotalAvailEnergy+TotalResEnergy
  # plot(days,E_habitat[,2]) # optional plots
  # plot(days,TotalAvailEnergy)
  # plot(days,TotalResEnergy)
  # plot(days,TotalEnergy)
  # plot(days,dailyDemand)
  # plot(days,TotalDemand)
  output<-data.frame(days,EnergyAvailable,EnergyReserve)
  write.csv(file="result.csv",output) # write data into file
  head(output) # show the top of the output file
  plot(days,energydeficit)
  # begin stacked plot
  AvailableE=matrix(t(E_food_type),nrow=(I*(totalsteps+1)),ncol=1)
  da1 <- data.frame(day=rep(1:(totalsteps+1),each=I),food=rep(foods),val=AvailableE)
  p<- ggplot(da1, aes(day,AvailableE))
  p+geom_area(aes(colour=food,fill=food),position='stack')
  # end of stacked plot

}

