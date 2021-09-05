#constructing network
#'@param Windows list,cells in each window
#'@param CD8TCellExp.trajectory expression profile for al cells
#'@param DynamicGene backgroud genes for constructing network
#'@param allTFs set the TFs
#'@param detectNum Threshold used to filter genes based on the number of expressed cells 
#'@param detectPro Threshold used to filter genes based on the proportion of expressed cells
#'@param meanExp Threshold used to filter genes based on the average expression 
DynNet_RF <- function(Windows, CD8TCellExp.trajectory, DynamicGene, allTFs, detectNum=10, detectPro=0.05, meanExp=1){
  #1. quality control for each window
  ExpofWindows <- lapply(Windows, function(x){
    windExp <- CD8TCellExp.trajectory[DynamicGene, x] 
    ind1 <- rowMeans(windExp>0) > detectPro 
    ind2 <- rowMeans(windExp) > meanExp 
    ind3 <- rowSums(windExp>0) > detectNum 
    windExp <- windExp[ind1&ind2&ind3, ]
    return(windExp)
  })
  
  #2. run GENIE3
  library(GENIE3)
  weightofWindows <- lapply(ExpofWindows, function(x){
    weightMat <- GENIE3(exprMatrix=x, regulators=intersect(allTFs, rownames(x)), 
                        targets = rownames(x), nCores=16)   
    weightdf  <- getLinkList(weightMat)                
    return(weightdf)
  })
  
  return(weightofWindows)
}


#Prune weak edges
#'@param Windows list,cells in each window
#'@param CD8TCellExp.trajectory expression profile for al cells
#'@param weightofWindows weights of all edges
#'@param nsd for each target, we filter the Connected edges based on the threshold, the weights mean of the edges plus n times standard deviation (mean + n * sd)
#'@param weightThr global filter. Weight is greater than the threshold (weightThr)
#'@param positivecor default 0, caculating the correlation between TF and target, if the correlation more than positivecor, we consider the edge is up regulated
DynNet_filter <- function(Windows, CD8TCellExp.trajectory, weightofWindows, 
                          nsd=1, weightThr, positivecor=0){
  #Step 1, Filter the TFs for each target gene
  netofWindows <- lapply(weightofWindows, function(x){
    targetweight <- split(x, as.character(x[,2])) 
    targetinteract <- lapply(targetweight, function(y){
      thr <- mean(y$weight)+(nsd*sd(y$weight)) 
      y[y$weight > thr, ]
    })
    t1 <- do.call("rbind", targetinteract)
    return(t1)
  })
  
  #Step 2, global filter. Weight is greater than the threshold weightThr
  netofWindows <- lapply(netofWindows, function(x){ x[x[,3]>weightThr,] })

  #Step 3, extract the up- and down- regulate information 
  netofWindows <- lapply(names(Windows), function(x){
    Exp <- CD8TCellExp.trajectory[, Windows[[x]]]
    spearmanCor <- apply(netofWindows[[x]], 1, function(y){
      cor(Exp[y[1], ], Exp[y[2], ], method = "spearman")
    })
    direct <- rep(0, length(spearmanCor))
    direct[spearmanCor > positivecor] <- 1
    direct[spearmanCor < (0-positivecor)] <- -1
    data.frame(netofWindows[[x]], spearmanCor, direct)
  })

  return(netofWindows)
}
