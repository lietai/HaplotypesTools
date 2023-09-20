#' Determine the block of SNPS that are considered in LD
#'
#' @param LinkageDisequilibriumMatrix matrix of SNPs x SNPs of linkage desiquilibrium
#' @param Threshold threshold to determine the block of snps that are considered in LD
#'
#' @return
#' @import igraph
#' @import data.table
#' @export
#'
#' @examples
#' library(HapToolkit)
#' HAP<-structure(list(
#' rs1 = c("C", "C", "C", "C", "T", "C", "C", "C", "C", "C", "C", "C"),
#' rs4 = c("G", "G", "G", "G", "A", "G", "G", "G", "G", "G", "G", "G"),
#' rs12 = c("G", "G", "G", "G", "A", "G", "G", "G", "G", "G", "G", "G"),
#' rs7 = c("A", "A", "A", "A", "C", "A", "A", "A", "A", "A", "A", "A"),
#' rs14 = c("G", "A", "A", "G", "A", "G", "A", "G", "A", "G", "G", "A"),
#' rs75 = c("C", "T", "T", "C", "T", "C", "T", "C", "T", "C", "C", "T"),
#' rs13 = c("C", "T", "T", "C", "T", "C", "T", "C", "T", "C", "C", "T"),
#' rs3 = c("A", "A", "A", "A", "A", "A", "G", "A", "A", "A", "A", "A"),
#' rs19 = c("T", "T", "T", "T", "T", "T", "C", "T", "T", "T", "T", "T"),
#' rs21 = c("T", "C", "C", "T", "T", "T", "T", "C", "C", "T", "C", "C")),
#' row.names = c("Ind1_A", "Ind1_B", "Ind2_A", "Ind2_B",
#'  "Ind3_A", "Ind3_B", "Ind4_A", "Ind4_B", "Ind5_A",
#'   "Ind5_B", "Ind6_A", "Ind6_B"), class = "data.frame")
#'
#' LD<-ComputeLinkagedisequilibriumMatrix(HAP,"Rsquare")
#' HaploBlocks<-DetermineBlocks(LD)
#'
DetermineBlocks<-function(LinkageDisequilibriumMatrix,Threshold=1){
  LD_DT<-data.table::data.table(LinkageDisequilibriumMatrix)
  SNP_ORDER<-colnames(LinkageDisequilibriumMatrix)
  LD_DT$Rs<-colnames(LinkageDisequilibriumMatrix)
  LD_MELTED<-data.table::melt(LD_DT,id="Rs")
  LD_MELTED<-data.frame(LD_MELTED)
  LD_MELTED$Rs<-factor(as.vector(LD_MELTED$Rs),levels=SNP_ORDER)
  LD_MELTED$variable<-factor(as.vector(LD_MELTED$variable),levels=SNP_ORDER)
  LD_MELTED$value[LD_MELTED$Rs==LD_MELTED$variable]<-1
  LD_GRAPH<-igraph::graph_from_data_frame(
    LD_MELTED[LD_MELTED$value>=Threshold,],
    directed = FALSE)
  BLOCKS<-igraph::decompose(LD_GRAPH)[order(unlist(lapply(igraph::decompose(LD_GRAPH),function(x){length(igraph::V(x))})))]
  BlocksSnps<-lapply(BLOCKS,function(x){igraph::V(x)$name})
  BlocksLengths<-unlist(lapply(BLOCKS,function(x){length(x)}))
  return(list(BlockSnp=BlocksSnps,BlocksLengths=BlocksLengths,BLOCKS=BLOCKS))
}


#' Make SNP Block Tree make a tree of haplotype block based on normalized mutual
#' information
#'
#' @template HAP
#' @import data.tree
#' @return
#' @export
#'
#' @examples
#' library(HapToolkit)
#' HAP<-structure(list(
#' rs1 = c("C", "C", "C", "C", "T", "C", "C", "C", "C", "C", "C", "C"),
#' rs4 = c("G", "G", "G", "G", "A", "G", "G", "G", "G", "G", "G", "G"),
#' rs12 = c("G", "G", "G", "G", "A", "G", "G", "G", "G", "G", "G", "G"),
#' rs7 = c("A", "A", "A", "A", "C", "A", "A", "A", "A", "A", "A", "A"),
#' rs14 = c("G", "A", "A", "G", "A", "G", "A", "G", "A", "G", "G", "A"),
#' rs75 = c("C", "T", "T", "C", "T", "C", "T", "C", "T", "C", "C", "T"),
#' rs13 = c("C", "T", "T", "C", "T", "C", "T", "C", "T", "C", "C", "T"),
#' rs3 = c("A", "A", "A", "A", "A", "A", "G", "A", "A", "A", "A", "A"),
#' rs19 = c("T", "T", "T", "T", "T", "T", "C", "T", "T", "T", "T", "T"),
#' rs21 = c("T", "C", "C", "T", "T", "T", "T", "C", "C", "T", "C", "C")),
#' row.names = c("Ind1_A", "Ind1_B", "Ind2_A", "Ind2_B",
#'  "Ind3_A", "Ind3_B", "Ind4_A", "Ind4_B", "Ind5_A",
#'   "Ind5_B", "Ind6_A", "Ind6_B"), class = "data.frame")
#'
#' HapTree<-MakeSNPBlockTree(HAP)
#' print(HapTree$DataTree)
#' library(igraph)
#' #library(networkD3)
#' #plot(as.igraph(HapTree$DataTree,directed=TRUE,direction="climb"))
#' #HapTreeNetwork <- ToDataFrameNetwork(HapTree$DataTree, "name")
#' #simpleNetwork(HapTreeNetwork[-3], fontSize = 12)
MakeSNPBlockTree<-function(HAP){
  AGGREGATED_SNPS<-list()
  Round<-1
  SNPS_ORIGINAL_ORDER<-colnames(HAP)
  SNP_ORDER<-c()
  LD_DF<-ComputeLinkagedisequilibriumMatrix(HAP,LinkageDesiquilibriumType = "NI")
  LD_DF<-data.frame(LD_DF)
  while(ncol(HAP)>1 & max(as.matrix(LD_DF))>0){
    HAPS<-igraph::decompose(igraph::graph_from_data_frame(which(as.matrix(LD_DF) == max(as.matrix(LD_DF)), arr.ind = TRUE),directed = FALSE))
    cat("R",Round,sum(unlist(lapply(HAPS,length))),max(as.matrix(LD_DF)),"\n")
    STEP_NAME<-colnames(HAP)
    for(CURRENT_HAP in 1:length(HAPS)){
      SNPS_TO_CONCAT<-STEP_NAME[sort(as.numeric(igraph::V(HAPS[[CURRENT_HAP]])$name))]
      SNP_ORDER<-c(SNP_ORDER,SNPS_TO_CONCAT)
      LocalHAPName<-paste("Hap_",Round,"_",CURRENT_HAP,sep='')
      HAP[,LocalHAPName]<-apply(HAP[,SNPS_TO_CONCAT],1,function(x){paste(x,collapse ="")})
      AGGREGATED_SNPS[[LocalHAPName]]<-SNPS_TO_CONCAT
      HAP[,SNPS_TO_CONCAT]<-NULL
      LD_DF[LocalHAPName,]<-0
      LD_DF[,LocalHAPName]<-0
      for(j in colnames(HAP)){
        if(LocalHAPName!=j){
          LD_DF[j,LocalHAPName]<-ComputeLinkagedisequilibrium(HAP,LocalHAPName,j,LinkageDesiquilibriumType ="NI")
        }
        LD_DF[,SNPS_TO_CONCAT]<-NULL
        LD_DF<-LD_DF[!rownames(LD_DF) %in% SNPS_TO_CONCAT,]
      }
    }
    Round<-Round+1
  }
  for(i in c(SNPS_ORIGINAL_ORDER,names(AGGREGATED_SNPS))){
    eval(parse(text=(paste(i,"<-data.tree::Node$new(\"",i,"\")"))))
  }

  for(i in names(AGGREGATED_SNPS)){
    for(j in AGGREGATED_SNPS[[i]]){
      eval(parse(text=(paste0(i,"$AddChildNode(",j,")"))))
    }
  }
  eval(parse(text=(paste0("Haplotypes_Dendo<-as.dendrogram(",LocalHAPName,")"))))
  return(list(DataTree=get(LocalHAPName),Dendo=get("Haplotypes_Dendo")))
}

