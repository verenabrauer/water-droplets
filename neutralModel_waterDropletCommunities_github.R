# R code belonging to the manuscript submitted  by Brauer et al. about
# "Selection in microbial islands creates a large core community with variable relative abundances"
# This codes simulates neutral drift and dispersal in the water droplet communities from the Pitch Lake.


# _______________________________________________
# definition:
# _______________________________________________
# ps: phyloseq object with one randomly selected droplet per sampling site; 
# includes an abundance table with amplicon sequence variants (ASVs) and a standardized 
# sample size of 20000 (i.e. number of individuals in a the water droplet community). 
# Rare ASVs occurring less then 4 times are removed from the abundance table.
# Abundance table with row=samples, columns=sequences.
require(phyloseq)


# _______________________________________________
# set parameters:
# _______________________________________________
# ps: must be provided
Size <- 2e4     # number of cells or sequences per droplet
NoDroplets <- nsamples(ps)   # number of droplets (24)
MaxGen <- 12*10000    # number of generations
GenerationTime <- 30 # Generation time in days; determines also the frequency of neutral events/replacements
Dispersal <- 0 # Fraction of neutral events due to dispersal (0 = complete dispersal limitation, 1 = full mixing)


# _______________________________________________
# sim_drift2()            
# _______________________________________________
# This function simulates neutral drift and dispersal with an infinite metacommunity

sim_drift2 <- function(pso, homogenize, size, nodroplets, maxgen, gentime, disp) {
  
  # pso         phyloseq object with absolute abundances
  # homogenize  logical; whether to create homogeneous starting communities or not
  # size        number of cells within droplet community
  # nodroplets  number of droplets to simulate
  # maxgen      number of generations to simulate
  # gentime     generation time
  # disp        fraction of neutral events due to dispersal (between 0 and 1; 0 = no dispersal, 1 = full mixing)
  
  
  # intergrated phyloseq-object with sum off all samples:
  nsamples <- nrow(otu_table(pso));   
  temp_data <- cbind(sample_data(pso),"group"=rep("sum", nsamples)) 
  sample_data(pso) <- temp_data                             
  pso.sum <- merge_samples(pso, "group"==group, fun = sum) 
  
  nasv <- ncol(otu_table(pso.sum))    # total number of ASVs accross all communities
  timesteps <- maxgen*gentime         # total number of time steps with neutral drift events
  stepsize <- size/gentime            # number of drift events per time step/ number of individuals to be replaced per generation
  if (stepsize < 1) {
    stepsize <- 1; timesteps <- maxgen*(size/stepsize)
  }
  metaS <- size*nodroplets            # size (number of individuals) of metacommunity
  
  
  # create metacommunity (pool (array) with  all sequences from all droplets):
  scom <- list()                      
  for (i in 1:nasv){
    scom[[i]] <-  rep(colnames(otu_table(pso.sum))[i],unname(otu_table(pso.sum)[,i]))
  }
  scom <- unlist(scom)               
  
  
  # randomly fill initial communities with sequences from metacommunity:
  droplet.composition <- character(size)
  #set.seed(seed)             # specify seed for exact reproducibility
  for (i in 1:size){
    droplet.composition[i] <- sample(scom, 1)
  }
  
  
  # determine whether initial communities differ or are exact replicates:
  if(homogenize==TRUE){        # equal droplet communities
    community.start <- list()
    for (k in 1:nodroplets){
      community.start[[k]] <- droplet.composition 
    }
  }
  else{                       # variable (real) droplet communities
    community.start <- list()
    for (k in 1:nodroplets){
      dcomp <- list()
      for (l in 1:nasv){
        dcomp[[l]] <-  rep(colnames(otu_table(pso))[l],unname(otu_table(pso)[k,l]))
      }
      dcomp <- unlist(dcomp)
      community.start[[k]] <- dcomp
      community.start[[k]] <- sample(community.start[[k]],size, replace=TRUE)
    }
  }
  
  
  
  # simulation of neutral drift and dispersal:
  community.evol <- community.start
  events <- sample(1:size, stepsize)               # individuals that undergo a neutral event
  lim <- round((1-disp)*length(events),0)
  
  # return number of time steps to skip until next dispersal event (for low dispersal rates):
  if (length(events)-lim>0) {gap<-0} else{gap <- round(1/(lim-round((1-disp)*length(events),3)),0)}; 
  
  if(timesteps>0){
    for (i in 1:timesteps) {
      if (disp == 0 | (gap>0 & i%%gap!=0)){        # dispersal = 0
        for (k in 1:nodroplets){
          community.evol[[k]][sample(1:size, stepsize)] <- community.evol[[k]][sample(1:size, stepsize)]
        }
      }
      else{                                        # dispersal > 0
        for (k in 1:nodroplets){
          events <- sample(1:size, stepsize)     
          lim <- round((1-disp)*length(events),0)
          if (gap>0 & i%%gap==0){                  # neutral drift and dispersal (low; one dispersal event after 'gap' times steps)
            intl <- events[1:(length(events)-1)]   # number of neutral events ocurring internally (within a droplet) = drift events
            extl <- events[length(events)]         # number of neutral events due to dispersal between droplets
            community.evol[[k]][intl] <- community.evol[[k]][sample(1:size, length(intl))]
            community.evol[[k]][extl] <- sample(scom,length(extl))
            rm(events,lim,intl,extl)
          }
          else{                                    # neutral drift and dispersal (every timestep):
            intl <- events[1:lim]                  # number of neutral events ocurring internally (within a droplet) = drift events
            extl <- events[(lim+1):length(events)] # number of neutral events due to dispersal between droplets
            community.evol[[k]][intl] <- community.evol[[k]][sample(1:size, length(intl))]
            community.evol[[k]][extl] <- sample(scom,length(extl))
            rm(events,lim,intl,extl)
          }
        }
      }
    }
    rm(timesteps,stepsize, metaS)
  }
  
  
  # reorganize data: create sequence table and combine with metadata
  list.start <- list(); list.evol <- list();
  for (k in 1:nodroplets){
    list.start[[k]] <- data.frame(table(community.start[[k]]))  # make a count table for each droplet (initial communities)
    list.start[[k]]$census <- "Start"                           # add a census column 
    list.start[[k]]$sample <- paste0("D",k)                     # add a sample name column
    ve.s <- as.character(list.start[[k]]$Var1)
    list.start[[k]]$asv <- unname(as.character(tax_table(pso.sum)[ve.s,7]))
    
    list.evol[[k]] <- data.frame(table(community.evol[[k]]))    # make a count table for each droplet (final communities)
    list.evol[[k]]$census <- "Evolved"                          # add a census column
    list.evol[[k]]$sample <- paste0("D",k)                      # add a sample name column
    ve.e <- as.character(list.evol[[k]]$Var1)
    list.evol[[k]]$asv <- unname(as.character(tax_table(pso.sum)[ve.e,7]))
    rm(ve.s,ve.e)
  }
  
  ## complete dataframe with initial and final communities:
  df.all <- rbind(list.start[[1]],list.evol[[1]])
  for (k in 2:nodroplets){
    df.all <- rbind(df.all, list.start[[k]], list.evol[[k]])
  }
  colnames(df.all)[1] <- "sequence"     # rename column 1
  colnames(df.all)[2] <- "abundance"    # rename column 2
  df.all$census_f <- factor(df.all$census, levels = c("Start","Evolved")) 
  
  return(df.all)
}



# _______________________________________________
# run neutral model
# _______________________________________________
rm(dfSim)

dfSim <- sim_drift2(pso=ps, homogenize=FALSE, size=Size, nodroplets=NoDroplets, 
                    maxgen=MaxGen, gentime=GenerationTime, disp=Dispersal)






