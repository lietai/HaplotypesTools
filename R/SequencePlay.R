
#' Compute hamming distance matrix between haplotype sequence
#'
#' @template HAP
#'
#' @import data.table
#' @import igraph
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
#' DIST<-ComputeHaplotypeSequenceDistanceGraph(HAP)
#' plot(DIST, layout=igraph::layout_with_kk(DIST,weights = igraph::E(DIST)$value))
#' #library(networkD3)
#' #simpleNetwork(as_long_data_frame(DIST)[,4:5])
ComputeHaplotypeSequenceDistanceGraph<-function(HAP){
  HaplotypesAllele<-table(apply(HAP,1,function(y){paste0(y,collapse="")}))
  TEST_LOCAL_LIST<-strsplit(names(HaplotypesAllele),'')
  names(TEST_LOCAL_LIST)<-names(HaplotypesAllele)
  HAP_HAMMING_DIST<-matrix(Inf,length(names(HaplotypesAllele)),length(names(HaplotypesAllele)))
  colnames(HAP_HAMMING_DIST)<-names(TEST_LOCAL_LIST)
  rownames(HAP_HAMMING_DIST)<-names(TEST_LOCAL_LIST)

  for(a in 1:length(names(TEST_LOCAL_LIST))){
    for(b in a:length(names(TEST_LOCAL_LIST))){
      i<-names(TEST_LOCAL_LIST)[a]
      j<-names(TEST_LOCAL_LIST)[b]
      #cat(i,j)
      HAP_HAMMING_DIST[i,j]<-sum(TEST_LOCAL_LIST[[i]]!=TEST_LOCAL_LIST[[j]])
    }
  }
  HAP_HAMMING_DIST<-data.frame(HAP_HAMMING_DIST)
  HAP_HAMMING_DIST$Hap_B<-rownames(HAP_HAMMING_DIST)
  HAP_HAMMING_DIST_DT<-data.table::melt(data.table::data.table(HAP_HAMMING_DIST),id="Hap_B")
  HAP_HAMMING_DIST_DT<-HAP_HAMMING_DIST_DT[order(HAP_HAMMING_DIST_DT$value),]
  HAP_HAMMING_DIST_DT<-HAP_HAMMING_DIST_DT[HAP_HAMMING_DIST_DT$Hap_B!=HAP_HAMMING_DIST_DT$variable,]
  HAP_HAMMING_DIST_DT$variable<-as.vector(HAP_HAMMING_DIST_DT$variable)

  HAP_SEQUENCE_GRAPH<-igraph::make_empty_graph(n=length(rownames(HAP_HAMMING_DIST)),directed=FALSE)
  igraph::V(HAP_SEQUENCE_GRAPH)$name<-rownames(HAP_HAMMING_DIST)
  while(igraph::components(HAP_SEQUENCE_GRAPH)$no>1 & dim(HAP_HAMMING_DIST_DT)[1]>0){
    HAP_SEQUENCE_GRAPH<-HAP_SEQUENCE_GRAPH+ igraph::edge(as.vector(unlist(HAP_HAMMING_DIST_DT[1,1])), as.vector(unlist(HAP_HAMMING_DIST_DT[1,2])), value = as.vector(unlist(HAP_HAMMING_DIST_DT[1,3])))
    if(!all(unlist(lapply(igraph::decompose(HAP_SEQUENCE_GRAPH),igraph::is_tree)))){
      HAP_SEQUENCE_GRAPH<-igraph::delete.edges(HAP_SEQUENCE_GRAPH,paste(as.vector(unlist(HAP_HAMMING_DIST_DT[1,1:2])),collapse ="|"))
    }
    HAP_HAMMING_DIST_DT<-HAP_HAMMING_DIST_DT[-1,]
  }
  return(HAP_SEQUENCE_GRAPH)
}


#' Title
#'
#' @param NSnps Number of SNPS to take account for this list
#' @param VCFStructure The structure obtained for this packages
#'
#' @return data frame for of the different combination of NSnps and the max
#' difference for all observed
#' @export
#'
#' @import data.table
#'
#'
#' @examples
#' library("HapToolkit")
#' VCF_file <- system.file("extdata",
#'  "TEST.vcf",
#'  package = "HapToolkit")
#'  FAM_file <- system.file("extdata",
#'  "TEST.fam",
#'  package = "HapToolkit")
#' I<-ReadPhasedVCF(VCF_file,FAM_file)
#' HAPS<-EnumerateAllHapsOfN(I,2);HAPS[order(HAPS$MaxDiff,decreasing = TRUE),]
#' start_time <- Sys.time()
#' HAP<-data.frame()
#' for (i in 2:7){
#'    cat(i,"\n")
#'   TEMP<-EnumerateAllHapsOfN(I,i)
#'   HAP<-rbind(HAP,TEMP)
#' }
#' end_time <- Sys.time()
#' end_time - start_time
#'# start_time <- Sys.time()
#'# HAP<-data.frame()
#'# ZeMaxDiff<-0
#'# for (i in 2:length(I$SNP_LIST)){
#'#   cat(i,"\n")
#'#   TEMP<-EnumerateAllHapsOfN(I,i)
#'#   if(max(TEMP$MaxDiff)<ZeMaxDiff){
#'#     break
#'#   }else{
#'#     ZeMaxDiff<-max(TEMP$MaxDiff)
#'#     HAP<-rbind(HAP,TEMP)
#'#   }
#'# }
#'# end_time <- Sys.time()
#'# end_time - start_time
EnumerateAllHapsOfN<-function(VCFStructure,NSnps=2){
  HAPLOTYPES<-data.frame(t(utils::combn(VCFStructure$SNP_LIST,NSnps)))
  HAPLOTYPES$MaxDiff<-apply(HAPLOTYPES,1,function(X){
    CASE_COUNTS<-table(
      apply(VCFStructure$CASES[,X],
            1,
            function(y){
              paste(y,sep="",collapse = "")
              }))
    CTRL_COUNTS<-table(
      apply(VCFStructure$CTRL[,X],
            1,
            function(y){
              paste(y,sep="",collapse = "")
              }))
    PRESENT_IN_CONTROLS<-names(CTRL_COUNTS)
    PRESENT_IN_CASES<-names(CASE_COUNTS)
    ABSENT_IN_CONTROLS<-PRESENT_IN_CASES[!PRESENT_IN_CASES %in% PRESENT_IN_CONTROLS]
    ABSENT_IN_CASES<-PRESENT_IN_CONTROLS[!PRESENT_IN_CONTROLS %in% PRESENT_IN_CASES]
    if(length(ABSENT_IN_CONTROLS)>=1){
      CTRL_COUNTS<-c(
        CTRL_COUNTS,
        stats::setNames(
          rep(0,length(ABSENT_IN_CONTROLS)),
          ABSENT_IN_CONTROLS))
    }
    if(length(ABSENT_IN_CASES)>=1){
      CASE_COUNTS<-c(
        CASE_COUNTS,
        stats::setNames(
          rep(0,length(ABSENT_IN_CASES)),
          ABSENT_IN_CASES))
    }
    max(abs(CASE_COUNTS-CTRL_COUNTS[names(CASE_COUNTS)]))
  })
  HAPLOTYPES$Combination<-apply(
    HAPLOTYPES[,colnames(HAPLOTYPES)[colnames(HAPLOTYPES)!="MaxDiff"]],1,
    function(x){paste(x,collapse="_")})
  return(HAPLOTYPES[,c("Combination","MaxDiff")])
}



