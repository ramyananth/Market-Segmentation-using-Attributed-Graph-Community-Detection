require(igraph)
require('proxy')

# Read values from command line
args <- commandArgs(trailingOnly = TRUE)
alpha = as.numeric(args[1])

# Read files
setwd("C:/Users/Ramya Ananth/Downloads/Attributed Graph Compunity Detection Project Resources-20190309")
attrib <- read.csv("/data/fb_caltech_small_attrlist.csv")
graph <- read.graph("/data/fb_caltech_small_edgelist.txt","edgelist")


# Computing Similarity matrix
SimilarityA <- as.matrix(simil(attrib,method = "cosine"))

# Computing Phase 1
initial<<-1:vcount(graph)
phase1 <- function(graph)
{
  iterations <- 15
  community <- 1:vcount(graph)
  temp_comm <- community
  initialMod <<-- Inf
  sim <- c()
  for(k in 1:nrow(SimilarityA))
  {
    sum <- 0
    for(m in 1:ncol(SimilarityA))
    {
      
      if(k!=m)
      {
        sum <- sum+SimilarityA[k,m]
      }
    }
    sim<-c(sim,sum)
  }
  continue <- TRUE
  while(continue)
  {
    temp_comm <- community
    for(i in 1:(length(community)-1))
    {
      max_delta <-- Inf
      max_j <- 0
      mod <- modularity(graph,community)
      for(j in (i+1):length(community))
      {
          if(i!=j)
          {
            new_community <- rep(community)
            new_community[i] <- j
            den <- length(unique(community))*length(community[community==1:length(community)]==TRUE)
            delta <- alpha*(modularity(graph,new_community)-mod)+((1-alpha)*sim[j])/den
            
            if(length(delta)!=0 && delta>max_delta &&delta>0)
            {
              max_delta <- delta
              max_j <- j
            }
          }
      }
      if(max_j!=0)
      {
        if(community[max_j]>community[i])
        {
          community[max_j] <- community[i]
        }
        else
        {
          community[i] <- community[max_j]
        }
      }

    }
    new_mod <- modularity(graph,community)
    iterations <- iterations-1
    if(initialMod<new_mod && iterations>0 && new_mod>0)
    {
      initialMod <<- new_mod
    }
    else
    {
      continue <- FALSE
      if(iterations>12 && alpha!=0)
      {
        community <- temp_comm
      }
    }
  }
  community
}

#Building Mapping to check for communities and building new graph with new vertices
mapping <- function(community)
{
  for(i in 1:length(initial))
  {
    if(is.na(community[initial[i]]))
    {
      t <- initial[i]
      print(c(initial[i],community[initial[i]],t))
    }
    initial[i] <<- community[initial[i]]
  }

  count <- 1
  mapping <- data.frame(matrix(0,nrow=length(unique(community)), ncol = 2))
  for(i in 1:length(community))
  {
    flag <- FALSE
    if(count>1)
    {
      for(j in 1:(count-1))
      {
        if(community[i]==mapping[j,1])
        {
          flag <- TRUE
          break
        }
      }
      if(flag==FALSE)
      {
        mapping[count,1] <- community[i]
        mapping[count,2] <- count
        count <- count+1
      }
    }
    else
    {
      mapping[1,1]<-community[i]
      mapping[1,2]<-count
      count<-count+1
    }
  }
  for(i in 1:length(community))
  {
    for(j in 1:nrow(mapping))
    {
      if(community[i]==mapping[j,1])
      {
        community[i]<-mapping[j,2]
        break
      }
    }
  }
  for(i in 1:length(initial))
  {
    for(j in 1:(count-1))
    {
      if(initial[i]==mapping[j,1])
      {
        initial[i]<<-mapping[j,2]
        break
      }
    }
  }
  community
}

# Updating Similarity Matrix
updatingSimilarityA<-function(community)
{
  for(i in 1:length(community))
  {
    temp<-SimilarityA[i,community[i]]+SimilarityA[community[i],i]
    SimilarityA[i,community[i]]<<-temp
    SimilarityA[community[i],i]<<-temp
  }
}

# Computing Phase 2
phase2 <- function()
{
  iterations_p2 <- 15

  while(iterations_p2>0 && vcount(graph)>1)
  {
    initialMod<<--Inf
    community<-phase1(graph)
    mapping<-mapping(community)
    print(length(unique(initial)))
    graph<-contract.vertices(graph,mapping)
    graph<-simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
    updatingSimilarityA(community)
    iterations_p2<-iterations_p2-1
  }
  
  # Write to the file
  
  fileName<-paste("./communities",alpha,sep="_")
  fileName<-paste(fileName,"txt",sep=".")
  fileptr<-file(fileName,"w")
  
  for(i in 1:max(initial))
  {
    ptr<-vector("numeric")
    for(j in 1:length(initial))
    {
      if(initial[j]==i){
        ptr<-append(ptr,j-1,after=length(ptr))
      }
    }
    cat(as.character(ptr),file=fileptr,sep = ",")
    cat("\n",file=fileptr)
  }
  close(fileptr)
}

phase2()