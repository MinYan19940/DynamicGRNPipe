###
#Example
###
# The example data can be downloaded at https://www.jianguoyun.com/p/DeYtp_AQ1bXiCRjXvYwE
#====0.input====
load(file = "clusterSig.RData")#genes used to construct cell trajectories
load(file = "GSE99254/CD8TCellExp.norm.RData") #expression profile and cell annotation

#====1. Construction of cell state transformation trajectory (slingshot)====
library(slingshot)
source("1.slingshot.R")
t.slingshot = slingshot_run(CD8TCellExp.norm,
                            clusterLabels=CellInfor$majorCluster, 
                            ordergene=unlist(clusterSig),
                            RMmethod="pca",
                            start.cluster="CD8_C1-LEF1") #
CellInfor.trajectory <- cbind.data.frame(CellInfor, t.slingshot$data) 
CD8TCellExp.trajectory <- CD8TCellExp.norm

#====2.slinding windows based on pseudotime and anotation of cells====
source("2.CellWindowing.R")
cells <- rownames(CellInfor.trajectory)[grep("1", CellInfor.trajectory$Branch)]  #cells in lineage1
cellInfor<- CellInfor.trajectory[cells, c("UniqueCell_ID", "majorCluster", "PseTime.curve1")]


#Calculate all intersections of pseudo-time density curves of cells in different states
C1_C2_cross <- densityintersection(a=cellInfor[cellInfor$majorCluster=="CD8_C1-LEF1","PseTime.curve1"], 
                                   b=cellInfor[cellInfor$majorCluster=="CD8_C2-CD28","PseTime.curve1"], 
                                   filename="C1_C2.pdf") 
C2_C4_cross <- densityintersection(a=cellInfor[cellInfor$majorCluster=="CD8_C2-CD28","PseTime.curve1"], 
                                   b=cellInfor[cellInfor$majorCluster=="CD8_C4-GZMK","PseTime.curve1"], 
                                   filename="C2_C4.pdf") 
C4_C5_cross <- densityintersection(a=cellInfor[cellInfor$majorCluster=="CD8_C4-GZMK","PseTime.curve1"], 
                                   b=cellInfor[cellInfor$majorCluster=="CD8_C5-ZNF683","PseTime.curve1"], 
                                   filename="C4_C5.pdf")
C5_C6_cross <- densityintersection(a=cellInfor[cellInfor$majorCluster=="CD8_C5-ZNF683","PseTime.curve1"], 
                                   b=cellInfor[cellInfor$majorCluster=="CD8_C6-LAYN","PseTime.curve1"], 
                                   filename="C5_C6.pdf")

#Manually select the intersection point according to the density plot
binpoint <- c(0, 16.4,26.29, 35.14, 43.51, max(cellInfor$PseTime.curve1))

#Two consecutive bins of cells as a window
#extract cells for each window
t1 <- cut(cellInfor$PseTime.curve1, binpoint) 
bincells <- split(cellInfor$UniqueCell_ID, t1)
W1 <- unlist(bincells[1:2])
W2 <- unlist(bincells[2:3])
W3 <- unlist(bincells[3:4])
W4 <- unlist(bincells[4:5])
Windows1 <- list(W1=W1, W2=W2, W3=W3, W4=W4)


#====3.constructing dynamic networks (GENIE3)======
source("3.constructionNetwork.R")
load("DynamicGene1.RData")
load("dorothea.RData")
#Construct a control network and calculate the control weight of each edge
weightofWindows <- DynNet_RF(Windows=Windows1, 
                             CD8TCellExp.trajectory=CD8TCellExp.trajectory, 
                             DynamicGene=DynamicGene1, #set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
                             allTFs=allTFs, #set regulators
                             detectNum=10, detectPro=0.05, meanExp=1 #Noise filtering threshold
                             )

#Filter the edges for each window
lineage1dynet <- DynNet_filter(Windows=Windows1, 
                               CD8TCellExp.trajectory=CD8TCellExp.trajectory, 
                               weightofWindows=weightofWindows, 
                               weightThr=0.02, 
                               nsd=2,  
                               positivecor=0, 
                               confidence=NULL)

#extract active edges for each window
Dynnet_active1 <- lapply(lineage1dynet, function(x){
  rownames(x) <- paste0(x[,1], "_", x[,2])
  x <- x[x[,"spearmanCor"]>0, ] 
  return(x)
})
names(Dynnet_active1) <- paste0("W", 1:length(Dynnet_active1))





