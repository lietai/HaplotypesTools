#' Title Compute Linkage disequilibrium between two snps or two locus
#'
#' @param data matrix
#'  ncols at least 2 cols
#'  n row = N individuals x 2 haplotypes
#' @param Snp1 SNP present int the data matrix
#' @param Snp2 SNP present int the data matrix
#' @param verbose print or not computation value False by default
#' @template LinkageDisequilibriumMethod
#'
#' @return single value of linkage disequilibrium
#' @export
#'
#' @examples
#'
#' Example<-structure(
#' c("G", "G", "A", "G", "G",
#'   "A", "G", "G", "A", "G",
#'   "C", "C", "C", "C", "C",
#'   "T", "C", "C", "C", "C",
#'   "G", "A", "A", "A", "A"),
#'   .Dim = c(5L, 5L),
#'   .Dimnames = list(
#'   c("Ind1_A","Ind1_B",
#'    "Ind2_A", "Ind2_B",
#'    "Ind3_A"),
#'    c("rs3", "rs1", "rs4", "rs44", "rs7" )))
#'ComputeLinkagedisequilibrium(Example,"rs3","rs44",LinkageDesiquilibriumType="Rsquare")
ComputeLinkagedisequilibrium<-function(
    data,
    Snp1,
    Snp2,
    LinkageDesiquilibriumType="D",
    verbose=FALSE){
    if(!LinkageDesiquilibriumType %in% c("D","Dprime","Rsquare","I","NI","NI2","NI3")){
      stop("LD unrecognized")
    }

    EFF<-table(data[,Snp1],data[,Snp2])
    F_EFF<-EFF/sum(EFF)
    F_A<-table(data[,Snp1])/sum(table(data[,Snp1]))
    F_B<-table(data[,Snp2])/sum(table(data[,Snp2]))
    F_A_Names<-names(F_A)
    F_B_Names<-names(F_B)
    I<-0
    for(A in F_A_Names){
      for(B in F_B_Names){
        if(verbose){cat(A,B,":",F_EFF[A,B],"\n")}
        H<-F_EFF[A,B]*log(F_EFF[A,B]/(F_A[A]*F_B[B]))
        if(!is.na(H)){
          I<-I+H
        }
      }
    }
    if(verbose){cat("Entropie:",I,"\n",sep=" ")}
    H_A<- sum(-F_A*log(F_A))
    H_B<- sum(-F_B*log(F_B))
    H_AB <- sum(-F_EFF*log(F_EFF),na.rm=TRUE)
    I<- H_A + H_B - H_AB
    NI <- (2*I)/(H_A + H_B)
    NI2 <- (I)/max(c(H_A,H_B))
    NI3 <- (I)/min(c(H_A,H_B))
    if(verbose){
      cat("Mutual Information:",I,"\n",sep=" ")
      cat("Normalized Mutual Information:",NI,"\n",sep=" ")
      cat("Normalized Mutual Information by max info:",NI2,"\n",sep=" ")
      cat("Normalized Mutual Information by min info:",NI3,"\n",sep=" ")
    }


    #if(sum(dim(EFF)==c(2,2))==2){
    i<-1
    j<-1
    EFF[i,j]/sum(EFF)-sum(EFF[,j])/sum(EFF)*sum(EFF[i,])/sum(EFF)

    i<-2
    j<-1
    -(EFF[i,j]/sum(EFF)-sum(EFF[,j])/sum(EFF)*sum(EFF[i,])/sum(EFF))
    i<-1
    j<-2
    -(EFF[i,j]/sum(EFF)-sum(EFF[,j])/sum(EFF)*sum(EFF[i,])/sum(EFF))
    i<-2
    j<-2
    D<-EFF[i,j]/sum(EFF)-sum(EFF[,j])/sum(EFF)*sum(EFF[i,])/sum(EFF)
    R_square<-(D*D)/(sum(EFF[,1])/sum(EFF)*sum(EFF[,2])/sum(EFF)*sum(EFF[1,])/sum(EFF)*sum(EFF[2,])/sum(EFF))

    if(D>0){
      Dmax<-min(sum(EFF[,1])/sum(EFF)*sum(EFF[2,])/sum(EFF),sum(EFF[,2])/sum(EFF)*sum(EFF[1,])/sum(EFF))
    }else{
      Dmax<-max(-(sum(EFF[,1])/sum(EFF)*sum(EFF[1,])/sum(EFF)),-(sum(EFF[,2])/sum(EFF)*sum(EFF[2,])/sum(EFF)))
    }
    if(verbose){
      cat("D:",D,"\n",sep=" ")
      cat("R square:",R_square,"\n",sep=" ")
      cat("D max:",Dmax,"\n",sep=" ")
    }


    Dprime<-D/Dmax
    if(LinkageDesiquilibriumType=="D"){
      return(D)
    }else if(LinkageDesiquilibriumType=="Dprime"){
      return(Dprime)
    }else if(LinkageDesiquilibriumType=="Rsquare"){
      return(R_square)
    }else if(LinkageDesiquilibriumType=="I"){
      return(I)
    }else if(LinkageDesiquilibriumType=="NI"){
      return(NI)
    }else if(LinkageDesiquilibriumType=="NI2"){
      return(NI2)
    }else if(LinkageDesiquilibriumType=="NI3"){
      return(NI3)
    }
    #}
    #else{return(-1)}
}

#' Compute Linkage disequilibrium matrix
#' from a haplotype matrix compute a SNP SNP LD matrix, the matrix is symetric
#'
#' @param HaplotypeMatrix matrix of haplotypes
#' @template LinkageDisequilibriumMethod
#' @return matrix SNPS x SNPS linkage disequilibrium
#' @export
#'
#' @examples
#'
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
#' rs21 = c("T", "C", "C", "T", "T", "T", "T", "C", "C", "T", "C", "C")
#' ),row.names = c("Ind1_A", "Ind1_B", "Ind2_A", "Ind2_B",
#'  "Ind3_A", "Ind3_B", "Ind4_A", "Ind4_B", "Ind5_A", "Ind5_B",
#'   "Ind6_A", "Ind6_B"), class = "data.frame")
#'
#' ComputeLinkagedisequilibriumMatrix(HAP,"Rsquare")
ComputeLinkagedisequilibriumMatrix<-function(HaplotypeMatrix,LinkageDesiquilibriumType){
  LD_MATRIX<-matrix(0,length(colnames(HaplotypeMatrix)),length(colnames(HaplotypeMatrix)))
  colnames(LD_MATRIX)<-colnames(HaplotypeMatrix)
  rownames(LD_MATRIX)<-colnames(HaplotypeMatrix)
  for(i in 1:length(colnames(HaplotypeMatrix))){
    for(j in i:length(colnames(HaplotypeMatrix))){
      #cat(i,"-",j,"\n")
      if(i!=j){
        LD_MATRIX[i,j]<-ComputeLinkagedisequilibrium(HaplotypeMatrix,i,j,LinkageDesiquilibriumType="NI")
        LD_MATRIX[j,i]<-ComputeLinkagedisequilibrium(HaplotypeMatrix,i,j,LinkageDesiquilibriumType="NI")
      }
    }
  }
  return(LD_MATRIX)
}


#' Title plot a Linkage desiquilibrium matrix as a ggplot2 graph
#'
#' @param LD_MATRIX matrix of linkage disequilibrium, suppose to be a square matrix
#' @template LinkageDisequilibriumMethod

#' @import data.table
#' @import ggplot2
#' @return G a ggplot2 graph of the linkage desiquilibrium
#' @export
#'
#' @examples
#' library(HapToolkit)
#' library(plotly)
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
#' ggplotly(PlotLinkagedisequilibriumMatrix(LD))
PlotLinkagedisequilibriumMatrix<-function(LD_MATRIX,LinkageDesiquilibriumType){
  LD_DF<-as.data.frame(LD_MATRIX)
  LD_DT<-data.table::data.table(LD_MATRIX)
  LD_DT$Rs<-colnames(LD_MATRIX)
  LD_MELTED<-data.table::melt(LD_DT,id="Rs")
  LD_MELTED<-data.frame(LD_MELTED)
  LD_MELTED$Rs<-factor(as.vector(LD_MELTED$Rs),levels=colnames(LD_MATRIX))
  LD_MELTED$variable<-factor(as.vector(LD_MELTED$variable),levels=colnames(LD_MATRIX))
  #https://stackoverflow.com/a/12429344/15246152
  #Rewrite your code to avoid non-standard evaluation. For ggplot2, this means using aes_string() instead of aes() (as described by Harlan)

  G<-ggplot2::ggplot(data=LD_MELTED,ggplot2::aes_string(x="Rs",y="variable",fill="value"))+
    ggplot2::geom_tile()+
    ggplot2::scale_y_discrete(position="right")+
    ggplot2::scale_fill_viridis_c(option="inferno",limits=c(0,1))
  return(G)
}



