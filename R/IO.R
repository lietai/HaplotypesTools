utils::globalVariables(c(":=", ".SD"))
#' ReadPhasedVCF
#'
#' Read a VCF and generate an haplotype matrix
#'
#' @param filename name of a phased vcf as for example
#' ##fileformat=VCFv4.2
#' #FILTER=<ID=PASS,Description="All filters passed">
#' #fileDate=24/05/2023 - 12:10:24
#' #source=shapeit4.1.3
#' #contig=<ID=1>
#' #INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#' #INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
#' #INFO=<ID=CM,Number=A,Type=Float,Description="Interpolated cM position">
#' #FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
#' #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  IndA       IndB   IndC
#' 1       944564  rs3128117       T       C       .       .       AC=269;AF=0.340506;CM=0.191843  GT      1|0     0|0     1|1
#' @param famfile .fam file used for Shape IT, as used in plink
#' @param chr chromosome of interest
#' @param start begening of the block of interest (position included)
#' @param end end of the block of interest
#' @import tool
#' @import data.table
#' @return list of data table
#' @export
#' @examples
#' library("HapToolkit")
#' VCF_file <- system.file("extdata",
#'  "TEST.vcf",
#'  package = "HapToolkit")
#'  FAM_file <- system.file("extdata",
#'  "TEST.fam",
#'  package = "HapToolkit")#'
#' I<-ReadPhasedVCF(VCF_file)
#' B<-ReadPhasedVCF(VCF_file,FAM_file)
#' C<-ReadPhasedVCF(VCF_file,FAM_file,chr=1,start=838555,end=891945)
#' filename<-VCF_file
#' famfile<-FAM_file
ReadPhasedVCF<-function(filename,famfile="",chr=0,start=0,end=0){
  verbose=TRUE
  ALT <- REF <- POS <- `#CHROM` <- NULL
  cat(filename,famfile)
  extension <- tools::file_ext(filename)
  if (extension != "vcf") {
    stop("File extension does not match '.vcf'")
  }
  if(famfile!=""){
    extension2 <- tools::file_ext(famfile)
    if (extension2 != "fam") {
      stop("File extension does not match '.fam'")
    }
  }
  VCF<-data.table::fread(filename)
  if(chr!=0 && start!=0 && end!=0){
    VCF<-VCF[`#CHROM`==chr][POS>=start][POS<=end]
  }
  if(verbose){ cat("data table done\n") }
  SAMPLES<-colnames(VCF)
  SAMPLES<-SAMPLES[-1:-9]
  VCF[,(SAMPLES ) := lapply(.SD,function(x){gsub(gsub(x,pattern = '1',  replacement =ALT),pattern = '0',replacement=REF)}), .SDcols = SAMPLES ,by=seq_len(nrow(VCF))]
  VCF<-data.frame(VCF)
  SAMPLES<-colnames(VCF)
  SAMPLES<-SAMPLES[-1:-9]
  if(verbose){ cat("Allele subtitued\n") }
  Haplotype=data.frame(
    rbind(
      as.matrix(
        data.frame(SAMPLES,"A")
        ),
      as.matrix(
        data.frame(SAMPLES,"B")
        )
      )
    )
  colnames(Haplotype)<-c("sample","strand")
  Haplotype$Haplotype<-""
  Haplotype$ID<-gsub("\\.","-",gsub("^X","",gsub("_.*","",Haplotype$sample)))
  for(i in SAMPLES){
    LEFT_NAME<-paste(i,"A",sep="_")
    RIGHT_NAME<-paste(i,"B",sep="_")
    VCF[,LEFT_NAME]<-gsub("\\|[ACGT]","",VCF[,i])
    VCF[,RIGHT_NAME]<-gsub("[ACGT]\\|","",VCF[,i])
    Haplotype[Haplotype$sample==i & Haplotype$strand=="A","Haplotype"]<-paste(VCF[,LEFT_NAME],collapse="")
    Haplotype[Haplotype$sample==i & Haplotype$strand=="B","Haplotype"]<-paste(VCF[,RIGHT_NAME],collapse="")
  }
  if(verbose){ cat("Haplotype computed\n") }
  SNP_LIST<-unlist(VCF$ID)
  H_VCF<-VCF
  H_VCF[,SAMPLES]<-NULL
  TH_VCF<-t(H_VCF)
  colnames(TH_VCF)<-TH_VCF["ID",]
  TH_VCF<-TH_VCF[-1:-9,]
  INFOS<-data.frame(ID=unlist(VCF$ID))
  rownames(INFOS)<-INFOS$ID
  con=file(filename,open="r")
  while (length(linn <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(length(grep("^##",linn,invert = TRUE))){
      break;
    }
    if(length(grep("^##INFO",linn))){
      TAG<-gsub(",.*$","",gsub("^.*ID=","",linn))
      cat(TAG,"\n")
      INFOS[,TAG]<-NA
    }

    #print(linn)
  }
  if(verbose){ cat("Get VCF infos\n") }
  close(con)
  INFOS_SPLIT<-strsplit(VCF$INFO,";")
  names(INFOS_SPLIT)<-VCF$ID
  for(SNP in names(INFOS_SPLIT)){
    lapply(INFOS_SPLIT[[SNP]],function(x){
      TAG<-gsub("=.*$","",x)
      VALUE<-gsub("^.*=","",x)
      INFOS[SNP,TAG]<<-VALUE
      })
  }
  if(verbose){ cat("Infos data frame computed\n") }
  INSTANCE<-list(
    VCF=data.table::data.table(VCF),
    Haplotype=Haplotype,
    SNP_LIST=SNP_LIST,
    HAP_SPLIT=TH_VCF,
    INFOS=INFOS)
  if(famfile!=""){
    cat(famfile)
    fam<-utils::read.table(famfile)
    fam$VCF_NAME<-paste(fam$V1,fam$V2,sep="_")
    fam$VCF_NAME<-gsub("-",".",fam$VCF_NAME)
    fam$A_NAME<-gsub("^([0-9])","X\\1",paste(fam$VCF_NAME,"A",sep="_"))
    fam$B_NAME<-gsub("^([0-9])","X\\1",paste(fam$VCF_NAME,"B",sep="_"))
    CasesHaP<-c(fam$A_NAME[fam$V6==2],fam$B_NAME[fam$V6==2])
    ControlsHaP<-c(fam$A_NAME[fam$V6==1],fam$B_NAME[fam$V6==1])
    CASES_HAP<-TH_VCF[CasesHaP,]
    CTRL_HAP<-TH_VCF[ControlsHaP,]
    INSTANCE$CASES<-CASES_HAP
    INSTANCE$CTRL<-CTRL_HAP
  }
  return(INSTANCE)
}


