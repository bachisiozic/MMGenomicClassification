##################################################################################
##################################################################################
##############################################################################1###
#                                                                                #
#        .     .     .                                               .     .     .
####     .     .     .   GISTIC - Identification of abnormal peaks   .     .     .     ####
#        .     .     .                                               .     .     .
#                                                                                #
##################################################################################
##################################################################################
##################################################################################

# # * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 
# #
# #                          Change the folder location (WD: working directory)
# #
# # * * * * * * * * * *         IMPORTANT         * * * * * * * * * * 
# WF="~/Google Drive/My Drive/Bachisio/Documents/Project/MM/Prediction_Model/GCP_MM/genomic/" # /where/is/the/folder/genomic
# 
# ##################################################################################
# ##################################################################################
# ##################################################################################
# ##################################################################################
# 
# 
# # Variables
# CYTOBAND.FILE=paste(WF,"/Reference/cytoBand.txt",sep="")
# ARMS=paste(WF,"/Reference/arm_hg19.txt",sep="")
# CHROM.LENGHT=paste(WF,"/Reference/hg19.chrom_sizes.txt",sep="")
# 
# GISTIC.FILE=paste(WF,"/Reference/Genomic_Coordinates_GISTICpeaks.txt",sep="")
# GISTIC.ARMS=paste(WF,"/Reference/arm_level.txt",sep="")
# 
# INPUT.PATH=paste(WF,"/input/cnv/",sep="")
# OUT.DIR=paste(WF,"/output/",sep="")
# PERCENTAGE=0.60 #(percentage of chromosome with gain/amp [for hyperdiploid chromosomes])


library(ggplot2)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library("GenomicRanges")
library(dplyr)
library(plyr)
library(ComplexHeatmap)
library(stringr)
library(stringi)


library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(annotate)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


library(IRanges)
library(data.table)
library(tidyr)
library(tibble)
library(MutationalPatterns)
library(deconstructSigs)
library(mmsig)
library(clonevol)

library(RColorBrewer)
library(pheatmap)
library(survival)
library(survminer)


#####################################################################################################################
#### . . . . . @ @ @ . . . . .       - FUNCTIONS for JCO Genomic Classification -       . . . . . @ @ @ . . . . .#### 
#####################################################################################################################


#### . . . -#### . . . - inputCNV()

#### . . . - Function to manage the CN file and create a consistent file with correct column names
####         DATA = 

inputCNV=function(DATA,
                  SAMPLE.COL,CHR.COL,START.COL,END.COL,TOTAL.COL,MINOR.COL){
  #### . . . - Changing column names in order to be consistent
  
  ############
  #### . . . - When the cancer fraction column is not defined 
  colnames(DATA)=ifelse(colnames(DATA)==SAMPLE.COL,"sample",
                        ifelse(colnames(DATA)==CHR.COL,"chrom",
                               ifelse(colnames(DATA)==START.COL,"start",
                                      ifelse(colnames(DATA)==END.COL,"end",
                                             ifelse(colnames(DATA)==TOTAL.COL,"total",
                                                    ifelse(colnames(DATA)==MINOR.COL,"minor",colnames(DATA)))))))
  
  #### . . . - Changing the variable type for start,end,total,minor,width,purity,ploidy,cancerfraction
  DATA$start=as.numeric(as.character(DATA$start))
  DATA$end=as.numeric(as.character(DATA$end))
  DATA$total=as.numeric(as.character(DATA$total))
  DATA$minor=as.numeric(as.character(DATA$minor))
  
  #### . . . - Formatting the chromosome column & keeping only the 24 chromosomes
  DATA$chrom=gsub("chr","",DATA$chrom)
  DATA$chrom=ifelse(DATA$chrom=="X",23,
                    ifelse(DATA$chrom=="Y",24,DATA$chrom))
  DATA=DATA[DATA$chrom %in% c(1:24),]
  
  #### . . . - Changing the variable type for chromosome column
  DATA$chrom=as.numeric(as.character(DATA$chrom))
  
  #### . . . - Sorting the object based on sample, chromosome, start
  DATA=DATA %>%
    arrange(sample,chrom,start)
  
  DATA=DATA[,c("sample","chrom","start","end","total","minor")]
  
  return(DATA)
}
  


#### . . . -#### . . . - FocalLargeMatrix()

#### . . . - Function to create the focal and large matrix for gain/amplification & biallelicloss/loss
####         CNV = 

FocalLargeMatrix=function(CNV,CYTOBAND.FILE,GISTIC.FILE,DATABASE,METADATA,FOCAL.LENGTH=5000000,OUT.DIR){
  ################################################################################################################################################################
  ######
  # . . . . .                                              Matrix with distinction focal and large
  ######
  ################################################################################################################################################################
  
  metadata=read.delim(METADATA,
                      stringsAsFactors = F)
  
  
  #####################
  ######## --- .                 Creation matrix with abnormal GISTIC peaks & CNV information
  #####################
  cytoband=read.delim(CYTOBAND.FILE,stringsAsFactors = F,header=F)
  colnames(cytoband)=c("X.chrom","chromStart","chromEnd","name","type")
  cytoband$cytoband=paste(gsub("chr","",cytoband$X.chrom),cytoband$name,sep="")
  
  gistic.peaks=read.delim(GISTIC.FILE,
                          stringsAsFactors = F)
  
  
  gistic.peaks.1=merge(gistic.peaks,cytoband,by="cytoband",all.x=T)
  gistic.peaks.1[gistic.peaks.1$type=="acen",]
  
  
  new.gistic=gistic.peaks
  new.gistic$gistic.seg.length=new.gistic$end-new.gistic$start
  
  # setwd(INPUT.PATH)
  # list.txt=list.files(patter="txt",full.names=T)
  # list.txt<-str_split_fixed(list.txt,"/",2)[,2]
  name.dataset=DATABASE
  
  all.cohort=NULL
  cnv=CNV
  cnv$type=name.dataset
  cnv$seg.length=cnv$end-cnv$start
  cnv$seg.type=ifelse(cnv$seg.length<=FOCAL.LENGTH,"focal","large")
  cnv$chrom=gsub("chr","",cnv$chrom)
  
  new.gistic$chr=gsub("chr","",new.gistic$chr)
    
  gr.cnv<- with(cnv, GRanges(chrom, IRanges(start=start, end=end)))
  values(gr.cnv) <- cnv[,c(1,5:ncol(cnv))]
    
  gr.gistic = with(new.gistic, GRanges(chr, IRanges(start=start, end=end)))
  values(gr.gistic) <- new.gistic[,c(4:ncol(new.gistic))]
    
  ranges <- merge(as.data.frame(gr.cnv),
                  as.data.frame(gr.gistic),
                  by="seqnames",
                  suffixes=c("A","B"))
  
  ranges_1=ranges[with(ranges, startB <= startA & endB <= endA & endB > startA),]
  ranges_2=ranges[with(ranges, startB >= startA & endB >= endA & startB < endA),]
  ranges_3=ranges[with(ranges, startB >= startA & endB <= endA),]
  ranges_4=ranges[with(ranges, startB < startA & endB > endA),]
  ranges.all<- rbind.data.frame(ranges_1, ranges_2, ranges_3, ranges_4)
    
  ranges.all=unique(ranges.all)
  ranges.all.1=ranges.all %>%
    arrange(sample,seqnames,startA)
    
  amp=ranges.all.1[ranges.all.1$lesion=="Amp",]
    # Gain Large 0.5
    # Gain Focal 1
    # Amp Large 1.5
    # Amp Focal 2
    # Normal => 0
  # amp$anomaly=ifelse(amp$total>=3.5 & amp$minor>=0.5 & amp$seg.type=="focal",2,
  #                    ifelse(amp$total>=3.5 & amp$minor>=0.5 & amp$seg.type=="large",1.5,
  #                           ifelse(amp$total>=2.5 & amp$total<3.5 & amp$minor>=0.5 & amp$seg.type=="focal",1,
  #                                  ifelse(amp$total>=2.5 & amp$total<3.5 & amp$minor>=0.5 & amp$seg.type=="large",0.5,0))))
  
  amp$anomaly=ifelse(amp$total>=3.5 & amp$seg.type=="focal",2,
                     ifelse(amp$total>=3.5 & amp$seg.type=="large",1.5,
                            ifelse(amp$total>=2.5 & amp$total<3.5 & amp$seg.type=="focal",1,
                                   ifelse(amp$total>=2.5 & amp$total<3.5 & amp$seg.type=="large",0.5,0))))
  
  
  amp=setDT(amp)
  amp.mtx=data.table::dcast(unique(amp[,c("sample","cytoband","anomaly")]),
                            formula=cytoband~sample,
                            toString,
                            value.var = "anomaly")
  # amp.mtx=data.table::dcast(as.data.table(unique(amp[,c("sample","cytoband","anomaly")])),
  #                           formula=cytoband~sample,
  #                           toString,
  #                           value.var = "anomaly")
  
  amp.mtx=as.data.frame(amp.mtx)
  rownames(amp.mtx)=amp.mtx$cytoband
  amp.mtx=amp.mtx[,-1]
  amp.mtx=as.matrix(amp.mtx)
  amp.mtx=ifelse(amp.mtx=="",NA,amp.mtx)
  amp.mtx=as.data.frame(amp.mtx)
  for(i in 1:length(colnames(amp.mtx))){
    for(z in 1:length(rownames(amp.mtx))){
      cel=amp.mtx[z,i]
      c1=str_split_fixed(cel,",",Inf)
      c1=as.matrix(c1)
      c1=gsub(" ","",c1)
      c1=as.data.frame(c1)
      for(zz in 1:length(colnames(c1))){
        c1[,zz]=as.numeric(c1[,zz])
      }
      c1$max=ifelse(is.na(c1),
                    NA,
                    max(c1,na.rm = T))
      amp.mtx[z,i]=c1$max
    }
  }
  missing.samples.amp=setdiff(unique(cnv$sample),colnames(amp.mtx))
  if(length(missing.samples.amp)>0){
    mis.amp.mtx=matrix(NA,nrow=length(rownames(amp.mtx)),ncol=length(missing.samples.amp))
    mis.amp.mtx=ifelse(is.na(mis.amp.mtx),0,mis.amp.mtx)
    colnames(mis.amp.mtx)=missing.samples.amp
    amp.mtx1=cbind(amp.mtx,mis.amp.mtx)
  }else{
    amp.mtx1=amp.mtx
  }
  rownames(amp.mtx1)=paste("Amp_",rownames(amp.mtx1),sep="")
  
  del=ranges.all.1[ranges.all.1$lesion=="Del",]
    # Biallelic Focal => 2
    # Biallelic Large => 1.5
    # Monoallelic Focal => 1
    # Monoallelic Large => 0.5
    # Normal => 0
  del$pr=ifelse(del$total>=2 & is.na(del$minor),0,"qqq")
  del$anomaly=ifelse(del$total<0.5 & del$seg.type=="focal",2,
                     ifelse(del$total<0.5 & del$seg.type=="large",1.5,
                            ifelse(del$total>=1.5 & del$minor<0.5 & del$seg.type=="focal",1,
                                   ifelse(del$total<1.5 & del$total>=0.5 & del$seg.type=="focal",1,
                                          ifelse(del$total>=1.5 & del$minor<0.5 & del$seg.type=="large",0.5,
                                                 ifelse(del$total<1.5 & del$total>=0.5 & del$seg.type=="large",0.5,0))))))
  
  del$anomaly=ifelse(is.na(del$anomaly),del$pr,del$anomaly)
    
  del=setDT(del)
  del.mtx=data.table::dcast(unique(del[,c("sample","cytoband","anomaly")]),
                            formula=cytoband~sample,
                            toString,
                            value.var = "anomaly")
  
  # del.mtx=data.table::dcast(as.data.table(unique(del[,c("sample","cytoband","anomaly")])),
  #                           formula=cytoband~sample,
  #                           toString,
  #                           value.var = "anomaly")
  del.mtx=as.data.frame(del.mtx)
  rownames(del.mtx)=del.mtx$cytoband
  del.mtx=del.mtx[,-1]
  del.mtx=as.matrix(del.mtx)
  del.mtx=ifelse(del.mtx=="",NA,del.mtx)
  del.mtx=as.data.frame(del.mtx)
  for(i in 1:length(colnames(del.mtx))){
    for(z in 1:length(rownames(del.mtx))){
      cel=del.mtx[z,i]
      c1=str_split_fixed(cel,",",Inf)
      c1=as.matrix(c1)
      c1=gsub(" ","",c1)
      c1=as.data.frame(c1)
      for(zz in 1:length(colnames(c1))){
        c1[,zz]=as.numeric(c1[,zz])
      }
      c1$max=ifelse(is.na(c1),
                    NA,
                    max(c1,na.rm = T))
      del.mtx[z,i]=c1$max
    }
  }
  missing.samples.del=setdiff(unique(cnv$sample),colnames(del.mtx))
  if(length(missing.samples.del)>0){
    mis.del.mtx=matrix(NA,nrow=length(rownames(del.mtx)),ncol=length(missing.samples.del))
    mis.del.mtx=ifelse(is.na(mis.del.mtx),0,mis.del.mtx)
    colnames(mis.del.mtx)=missing.samples.del
    del.mtx1=cbind(del.mtx,mis.del.mtx)
  }else{
    del.mtx1=del.mtx
  }
  rownames(del.mtx1)=paste("Del_",rownames(del.mtx1),sep="")
    
  all=rbind(amp.mtx1,del.mtx1)
  t.all=t(all)
  t.all=as.data.frame(t.all)
  for(jj in 1:length(colnames(t.all))){
    t.all[,jj]=as.numeric(as.character(t.all[,jj]))
  }
    
  mis.peak=setdiff(new.gistic$code,colnames(t.all))
  if(length(mis.peak)!=0){
    for(jjj in length(mis.peak)){
      t.all[,paste(mis.peak[jjj],"",sep="")]=0
    }
  }else{
    t.all=t.all
  }
  all.cohort=rbind(all.cohort,t.all)

  #### . . . when some peaks are missing
  pr=setdiff(new.gistic$code,unique(ranges.all$code))
  if(length(pr)>=1){
    pr.m=matrix(NA,nrow=length(rownames(all.cohort)),ncol=length(pr))
    colnames(pr.m)=pr
    rownames(pr.m)=rownames(all.cohort)
    all.cohort=cbind(all.cohort,pr.m)
  }
  
  all.cohort$database=DATABASE
  
  write.table(all.cohort,
              paste(OUT.DIR,"/CNV_",DATABASE,"matrix_focal_large_5Mb.txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
}



#### . . . - SegmentCollapsing()

#### . . . - Function to collapse the consecutive segments with the same total and minor 
####         copy_number = CN file

SegmentCollapsing=function(copy_number,DATABASE){
  new.DATA=NULL
  cnv=copy_number
  samples=unique(cnv$sample)
  for (i in 1:length(samples))
  {
    cnv.s=cnv[cnv$sample==samples[i],]
    chroms=unique(cnv.s$chrom)
    cnv.s$total=round(cnv.s$total,0)
    cnv.s$minor=round(cnv.s$minor,0)
    for(z in 1:length(chroms)){
      chr<-cnv.s[cnv.s$chrom==chroms[z],]
      chr$code=paste(chr$total,chr$minor,sep="_")
      
      a<-rle(chr$code)$lengths
      chr$first<-unlist(lapply(a, function(x) seq(1,x)))
      chr$last<-unlist(lapply(a, function (x) seq(x,1,-1)))
      
      chr1<-chr[chr$first==1 | chr$last==1,]
      
      a1<-rle(chr1$code)$lengths
      chr1$first2<-unlist(lapply(a1, function(x) seq(1,x)))
      chr1$last2<-unlist(lapply(a1, function (x) seq(x,1,-1)))
      
      ##_____ 
      chr1_start<-subset(chr1,chr1$first2==1)
      head(chr1_start)
      chr1_end<-subset(chr1,chr1$last2==1)
      head(chr1_end)
      colnames(chr1_end)<-c("sample","chrom1","start1","end1","total1","minor1","code","first1","last1","first21","last21")
      head(chr1_end)
      dim(chr1_start)
      dim(chr1_end)
      ##_____ 
      df1<-data.frame(chr1_start$sample,chr1_start$chrom,chr1_start$start,chr1_end$end1,chr1_end$total1,chr1_end$minor1,
                      chr1_end$code,chr1_end$first1)
      colnames(df1)<-c("sample","chrom","start","end","total","minor","code","Num_Segments")
      df1$type=DATABASE
      new.DATA=rbind(new.DATA,df1)
    }
  }
  return(new.DATA)
}



#### . . . - HYidentification()

#### . . . - Function to identify the hyperdiploid patients
####         CNV = 

HYidentification=function(CNV,CYTOBAND.FILE,ARMS,CHROM.LENGHT,GISTIC.ARMS,DATABASE,OUT.DIR,PERCENTAGE=0.6){
  ##################################################################################1
  #     .     .     .                                                   .     .     .
  #     .     .     .   PredictionModel - Identification hyperdiploid   .     .     .
  #     .     .     .                                                   .     .     .
  #B###############################################################################Z#
  
  
  ################################################################################################################################################################
  ######
  # . . . . .                                              ARMS
  ######
  ################################################################################################################################################################
  
  
  
  #####################
  ######## --- .                 Creation matrix with abnormal GISTIC arms & CNV information    
  #####################
  #####################
  ######## --- .                 CoMMpass, Bolli, MRCxi, MSK and UAMS   
  #####################
  cytoband=read.delim(CYTOBAND.FILE,stringsAsFactors = F,header=F)
  colnames(cytoband)=c("X.chrom","chromStart","chromEnd","name","type")
  cytoband$cytoband=paste(gsub("chr","",cytoband$X.chrom),cytoband$name,sep="")
  gistic.arms=read.delim(GISTIC.ARMS,
                         stringsAsFactors = F)
  
  new.gistic=gistic.arms
  
  name.dataset=DATABASE
  
  #### . . . - Segment collapsing with the same tot and min (based on code variable)   
  cnv.n=SegmentCollapsing(CNV,DATABASE = name.dataset)
  
  cnv.gain=cnv.n[cnv.n$total>=2.5 & cnv.n$minor>=0.5,]
  cnv.gain$seg.length=cnv.gain$end-cnv.gain$start
  cnv.gain1=cnv.gain
  
  cnv.gain=cnv.gain[complete.cases(cnv.gain$chrom),]
  
  #### . . . - Add arms information
  arms=read.delim(ARMS,stringsAsFactors = F)
  arms$arm.length=arms$End_Position-arms$Start_Position
  
  #### . . . - Creating the dataset with all information
  gr.cnv<- with(cnv.gain, GRanges(chrom, IRanges(start=start, end=end)))
  values(gr.cnv) <- cnv.gain[,c(1,5:ncol(cnv.gain))]
  
  gr.arms = with(arms, GRanges(Chromosome, IRanges(start=Start_Position, end=End_Position)))
  values(gr.arms) <- arms[,c(1,5)]
  
  ranges <- merge(as.data.frame(gr.cnv),
                  as.data.frame(gr.arms),
                  by="seqnames",
                  suffixes=c("A","B"))
  
  ranges_1=ranges[with(ranges, startB <= startA & endB <= endA & endB > startA),]
  ranges_2=ranges[with(ranges, startB >= startA & endB >= endA & startB < endA),]
  ranges_3=ranges[with(ranges, startB >= startA & endB <= endA),]
  ranges_4=ranges[with(ranges, startB < startA & endB > endA),]
  ranges.all<- rbind.data.frame(ranges_1, ranges_2, ranges_3, ranges_4)
  
  ranges.all=unique(ranges.all)
  ranges.all.1=ranges.all %>%
    arrange(type,sample,seqnames,startA)
  head(ranges.all.1)
  
  cnv.gain=ranges.all.1[,c("sample","seqnames","startA","endA","total","minor","code","Num_Segments","type","startB","endB","Arm","arm.length")]
  colnames(cnv.gain)=c(c("sample","chr","start","end","tot","min","code","Num_Segments","type","start.arm","end.arm","arm","arm.length"))
  
  #cnv.gain$p.q=substring(cnv.gain$arm,2) # keep the second character
  cnv.gain$p.q=stri_sub(cnv.gain$arm,-1)  # keep the last character
  
  cnv.gain$seg.length=ifelse(cnv.gain$end <= cnv.gain$end.arm & cnv.gain$start >= cnv.gain$start.arm,
                             cnv.gain$end-cnv.gain$start,
                             ifelse(cnv.gain$end > cnv.gain$end.arm & cnv.gain$start >= cnv.gain$start.arm,
                                    cnv.gain$end.arm-cnv.gain$start,
                                    ifelse(cnv.gain$start < cnv.gain$start.arm,
                                           cnv.gain$end-cnv.gain$start.arm,NA)))
  
  #### . . . - Chromosome identification with more than CN 2 in >60% of the chromosome (and no.hyperdiploid)
  chrom.length=read.delim(CHROM.LENGHT,
                          stringsAsFactors = F)
  f=c(23,155270560)
  m=c(24,59373566)
  chrom.length=rbind(chrom.length,f,m)
  
  cnv.gain$chr=as.character(cnv.gain$chr)
  
  hyperdiploid=NULL
  no.hyperdiploid=NULL
  samples=unique(cnv.n$sample)
  for(ii in 1:length(samples)){
    new.s=cnv.gain[cnv.gain$sample==samples[ii],]
    chroms=unique(new.s$chr)
    ARMS=unique(new.s$arm)
    if(length(new.s$sample)>0){
      for(zz in 1:length(chroms)){
        chr<-new.s[new.s$chr==chroms[zz],]
        for(zzz in 1:length(ARMS)){
          arm<-new.s[new.s$arm==ARMS[zzz],]
          sam=samples[ii]
          cro=unique(arm$chr)
          CN=median(as.numeric(as.character(arm$tot)),na.rm = T)
          max.cn=max(as.numeric(as.character(arm$tot)),na.rm = T)
          min.cn=min(as.numeric(as.character(arm$tot)),na.rm = T)  
          type=unique(arm$type)
          length.seg=sum(as.numeric(as.character(arm$seg.length)))
          p.q=unique(arm$p.q)
          seg.gain.prop=as.numeric(as.character(length.seg))/as.numeric(as.character((arms[arms$Arm==ARMS[zzz],"arm.length"])))
          line=c(sam,cro,CN,max.cn,min.cn,type,length.seg,p.q,seg.gain.prop)
          hyperdiploid=rbind(hyperdiploid,line)
        }
      }
    }else{
      no.hyperdiploid=rbind(no.hyperdiploid,samples[ii])
    }
  }
  
  hyper=as.data.frame(hyperdiploid)
  
  rownames(hyper)=c(1:length(hyper$V1))
  colnames(hyper)=c("sample","chromosome","CN.median","CN.max","CN.min","type","length.seg","arm","seg.gain.prop")
  
  hyper$seg.gain.prop=as.numeric(hyper$seg.gain.prop)
  hyper=unique(hyper)
  
  samples=unique(hyper$sample)
  hyper.new=NULL
  for(ii in 1:length(samples)){
    new.s=hyper[hyper$sample==samples[ii],]
    chroms=unique(new.s$chr)
    for(zz in 1:length(chroms)){
      chr<-new.s[new.s$chr==chroms[zz],]
      sam=samples[ii]
      cro=unique(chr$chr)
      CN.p=chr[chr$arm=="p","CN.median"]
      CN.p[length(CN.p)==0]=NA
      max.cn.p=chr[chr$arm=="p","CN.max"]
      max.cn.p[length(max.cn.p)==0]=NA
      min.cn.p=chr[chr$arm=="p","CN.min"]
      min.cn.p[length(min.cn.p)==0]=NA
      length.seg.p=chr[chr$arm=="p","length.seg"]
      length.seg.p[length(length.seg.p)==0]=NA
      arm.p=chr[chr$arm=="p","arm"]
      arm.p[length(arm.p)==0]=NA
      seg.gain.prop.p=chr[chr$arm=="p","seg.gain.prop"]
      seg.gain.prop.p[length(seg.gain.prop.p)==0]=NA
      
      CN.q=chr[chr$arm=="q","CN.median"]
      CN.q[length(CN.q)==0]=NA
      max.cn.q=chr[chr$arm=="q","CN.max"]
      max.cn.q[length(max.cn.q)==0]=NA
      min.cn.q=chr[chr$arm=="q","CN.min"]
      min.cn.q[length(min.cn.q)==0]=NA
      length.seg.q=chr[chr$arm=="q","length.seg"]
      length.seg.q[length(length.seg.q)==0]=NA
      arm.q=chr[chr$arm=="q","arm"]
      arm.q[length(arm.q)==0]=NA
      seg.gain.prop.q=chr[chr$arm=="q","seg.gain.prop"]
      seg.gain.prop.q[length(seg.gain.prop.q)==0]=NA
      
      type=unique(chr$type)
      gain.length=sum(as.numeric(as.character(chr$length.seg)))
      perc.gain.chr=as.numeric(as.character(gain.length))/(as.numeric(as.character(chrom.length[chrom.length$Offset==chroms[zz],"Size"])))
      line.p=c(sam,cro,CN.p,max.cn.p,min.cn.p,type,length.seg.p,arm.p,seg.gain.prop.p,gain.length,perc.gain.chr)
      line.q=c(sam,cro,CN.q,max.cn.q,min.cn.q,type,length.seg.q,arm.q,seg.gain.prop.q,gain.length,perc.gain.chr)
      hyper.new=rbind(hyper.new,line.p,line.q)
    }
  }
  
  
  hyper.new.1=as.data.frame(hyper.new)
  rownames(hyper.new.1)=c(1:length(hyper.new.1$V1))
  colnames(hyper.new.1)=c("sample","chromosome","CN.median","CN.max","CN.min","type","length.seg","arm","seg.gain.prop","gain.length","perc.gain.chr")
  
  hyper.new.2=hyper.new.1[complete.cases(hyper.new.1$CN.median),]
  hyper.new.2$arm.gain=ifelse(hyper.new.2$seg.gain.prop>PERCENTAGE,1,0)
  hyper.new.2$chr.gain=ifelse(hyper.new.2$perc.gain.chr>PERCENTAGE,1,0)
  head(hyper.new.2)
  
  hyper=hyper.new.2[hyper.new.2$chr.gain==1,]
  
  #### . . . - Hyperdiploid and DoubleGenome identification  (and no.hyperdiploid)
  chroms.hyper=c(3,5,7,9,11,15,18,19,21)
  samples=unique(hyper$sample)
  hyper.1=NULL
  no.hyperdiploid.1=NULL
  for(iii in 1:length(samples)){
    hyper.s=hyper[hyper$sample==samples[iii],]
    if(sum(hyper.s[,"chr.gain"])>20){
      sam=unique(hyper.s$sample)
      type=unique(hyper.s$type)
      whole.median.CN=median(as.numeric(as.character(hyper.s$CN.median)),na.rm = T)
      
      # . p arm
      chr3.p=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="p",]
      seg.gain.pro.chr3.p=ifelse(chr3.p$chromosome==3,chr3.p$seg.gain.prop,NA)
      seg.gain.pro.chr3.p[length(seg.gain.pro.chr3.p)==0]=NA
      chr5.p=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="p",]
      seg.gain.pro.chr5.p=ifelse(chr5.p$chromosome==5,chr5.p$seg.gain.prop,NA)
      seg.gain.pro.chr5.p[length(seg.gain.pro.chr5.p)==0]=NA
      chr7.p=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="p",]
      seg.gain.pro.chr7.p=ifelse(chr7.p$chromosome==7,chr7.p$seg.gain.prop,NA)
      seg.gain.pro.chr7.p[length(seg.gain.pro.chr7.p)==0]=NA
      chr9.p=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="p",]
      seg.gain.pro.chr9.p=ifelse(chr9.p$chromosome==9,chr9.p$seg.gain.prop,NA)
      seg.gain.pro.chr9.p[length(seg.gain.pro.chr9.p)==0]=NA
      chr11.p=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="p",]
      seg.gain.pro.chr11.p=ifelse(chr11.p$chromosome==11,chr11.p$seg.gain.prop,NA)
      seg.gain.pro.chr11.p[length(seg.gain.pro.chr11.p)==0]=NA
      chr15.p=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="p",]
      seg.gain.pro.chr15.p=ifelse(chr15.p$chromosome==15,chr15.p$seg.gain.prop,NA)
      seg.gain.pro.chr15.p[length(seg.gain.pro.chr15.p)==0]=NA
      
      chr18.p=hyper.s[hyper.s$chromosome==18 & hyper.s$arm=="p",]
      seg.gain.pro.chr18.p=ifelse(chr18.p$chromosome==18,chr18.p$seg.gain.prop,NA)
      seg.gain.pro.chr18.p[length(seg.gain.pro.chr18.p)==0]=NA
      
      chr19.p=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="p",]
      seg.gain.pro.chr19.p=ifelse(chr19.p$chromosome==19,chr19.p$seg.gain.prop,NA)
      seg.gain.pro.chr19.p[length(seg.gain.pro.chr19.p)==0]=NA
      chr21.p=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="p",]
      seg.gain.pro.chr21.p=ifelse(chr21.p$chromosome==21,chr21.p$seg.gain.prop,NA)
      seg.gain.pro.chr21.p[length(seg.gain.pro.chr21.p)==0]=NA
      # . q arm
      chr3.q=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="q",]
      seg.gain.pro.chr3.q=ifelse(chr3.q$chromosome==3,chr3.q$seg.gain.prop,NA)
      seg.gain.pro.chr3.q[length(seg.gain.pro.chr3.q)==0]=NA
      chr5.q=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="q",]
      seg.gain.pro.chr5.q=ifelse(chr5.q$chromosome==5,chr5.q$seg.gain.prop,NA)
      seg.gain.pro.chr5.q[length(seg.gain.pro.chr5.q)==0]=NA
      chr7.q=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="q",]
      seg.gain.pro.chr7.q=ifelse(chr7.q$chromosome==7,chr7.q$seg.gain.prop,NA)
      seg.gain.pro.chr7.q[length(seg.gain.pro.chr7.q)==0]=NA
      chr9.q=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="q",]
      seg.gain.pro.chr9.q=ifelse(chr9.q$chromosome==9,chr9.q$seg.gain.prop,NA)
      seg.gain.pro.chr9.q[length(seg.gain.pro.chr9.q)==0]=NA
      chr11.q=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="q",]
      seg.gain.pro.chr11.q=ifelse(chr11.q$chromosome==11,chr11.q$seg.gain.prop,NA)
      seg.gain.pro.chr11.q[length(seg.gain.pro.chr11.q)==0]=NA
      chr15.q=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="q",]
      seg.gain.pro.chr15.q=ifelse(chr15.q$chromosome==15,chr15.q$seg.gain.prop,NA)
      seg.gain.pro.chr15.q[length(seg.gain.pro.chr15.q)==0]=NA
      
      chr18.q=hyper.s[hyper.s$chromosome==18 & hyper.s$arm=="q",]
      seg.gain.pro.chr18.q=ifelse(chr18.q$chromosome==18,chr18.q$seg.gain.prop,NA)
      seg.gain.pro.chr18.q[length(seg.gain.pro.chr18.q)==0]=NA
      
      chr19.q=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="q",]
      seg.gain.pro.chr19.q=ifelse(chr19.q$chromosome==19,chr19.q$seg.gain.prop,NA)
      seg.gain.pro.chr19.q[length(seg.gain.pro.chr19.q)==0]=NA
      chr21.q=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="q",]
      seg.gain.pro.chr21.q=ifelse(chr21.q$chromosome==21,chr21.q$seg.gain.prop,NA)
      seg.gain.pro.chr21.q[length(seg.gain.pro.chr21.q)==0]=NA
      
      perc.gain.chr3=unique(hyper.s[hyper.s$chromosome==3,"perc.gain.chr"])
      perc.gain.chr3[length(perc.gain.chr3)==0]=NA
      perc.gain.chr5=unique(hyper.s[hyper.s$chromosome==5,"perc.gain.chr"])
      perc.gain.chr5[length(perc.gain.chr5)==0]=NA
      perc.gain.chr7=unique(hyper.s[hyper.s$chromosome==7,"perc.gain.chr"])
      perc.gain.chr7[length(perc.gain.chr7)==0]=NA
      perc.gain.chr9=unique(hyper.s[hyper.s$chromosome==9,"perc.gain.chr"])
      perc.gain.chr9[length(perc.gain.chr9)==0]=NA
      perc.gain.chr11=unique(hyper.s[hyper.s$chromosome==11,"perc.gain.chr"])
      perc.gain.chr11[length(perc.gain.chr11)==0]=NA
      perc.gain.chr15=unique(hyper.s[hyper.s$chromosome==15,"perc.gain.chr"])
      perc.gain.chr15[length(perc.gain.chr15)==0]=NA
      
      perc.gain.chr18=unique(hyper.s[hyper.s$chromosome==18,"perc.gain.chr"])
      perc.gain.chr18[length(perc.gain.chr18)==0]=NA
      
      perc.gain.chr19=unique(hyper.s[hyper.s$chromosome==19,"perc.gain.chr"])
      perc.gain.chr19[length(perc.gain.chr19)==0]=NA
      perc.gain.chr21=unique(hyper.s[hyper.s$chromosome==21,"perc.gain.chr"])
      perc.gain.chr21[length(perc.gain.chr21)==0]=NA
      
      
      double.chrs=sum(hyper.s[,"arm.gain"]) # Number of gained/amplified arms in the sample
      Cytogenetic="Double_Genomes"
      line.double.p=c(sam,type,whole.median.CN,seg.gain.pro.chr3.p,seg.gain.pro.chr5.p,seg.gain.pro.chr7.p,seg.gain.pro.chr9.p,
                      seg.gain.pro.chr11.p,seg.gain.pro.chr15.p,
                      seg.gain.pro.chr18.p,
                      seg.gain.pro.chr19.p,seg.gain.pro.chr21.p,double.chrs,Cytogenetic,"p",
                      perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,
                      perc.gain.chr18,
                      perc.gain.chr19,perc.gain.chr21)
      line.double.q=c(sam,type,whole.median.CN,seg.gain.pro.chr3.q,seg.gain.pro.chr5.q,seg.gain.pro.chr7.q,seg.gain.pro.chr9.q,
                      seg.gain.pro.chr11.q,seg.gain.pro.chr15.q,
                      seg.gain.pro.chr18.q,
                      seg.gain.pro.chr19.q,seg.gain.pro.chr21.q,double.chrs,Cytogenetic,"q",
                      perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,
                      perc.gain.chr18,
                      perc.gain.chr19,perc.gain.chr21)
      hyper.1=rbind(hyper.1,line.double.p,line.double.q)
    }else if(sum(hyper.s[hyper.s$chromosome %in% chroms.hyper,"chr.gain"])>=1){
      sam=unique(hyper.s$sample)
      type=unique(hyper.s$type)
      whole.median.CN=median(as.numeric(as.character(hyper.s$CN.median)),na.rm = T)
      
      # . p arm
      chr3.p=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="p",]
      seg.gain.pro.chr3.p=ifelse(chr3.p$chromosome==3,chr3.p$seg.gain.prop,NA)
      seg.gain.pro.chr3.p[length(seg.gain.pro.chr3.p)==0]=NA
      chr5.p=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="p",]
      seg.gain.pro.chr5.p=ifelse(chr5.p$chromosome==5,chr5.p$seg.gain.prop,NA)
      seg.gain.pro.chr5.p[length(seg.gain.pro.chr5.p)==0]=NA
      chr7.p=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="p",]
      seg.gain.pro.chr7.p=ifelse(chr7.p$chromosome==7,chr7.p$seg.gain.prop,NA)
      seg.gain.pro.chr7.p[length(seg.gain.pro.chr7.p)==0]=NA
      chr9.p=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="p",]
      seg.gain.pro.chr9.p=ifelse(chr9.p$chromosome==9,chr9.p$seg.gain.prop,NA)
      seg.gain.pro.chr9.p[length(seg.gain.pro.chr9.p)==0]=NA
      chr11.p=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="p",]
      seg.gain.pro.chr11.p=ifelse(chr11.p$chromosome==11,chr11.p$seg.gain.prop,NA)
      seg.gain.pro.chr11.p[length(seg.gain.pro.chr11.p)==0]=NA
      chr15.p=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="p",]
      seg.gain.pro.chr15.p=ifelse(chr15.p$chromosome==15,chr15.p$seg.gain.prop,NA)
      seg.gain.pro.chr15.p[length(seg.gain.pro.chr15.p)==0]=NA
      
      chr18.p=hyper.s[hyper.s$chromosome==18 & hyper.s$arm=="p",]
      seg.gain.pro.chr18.p=ifelse(chr18.p$chromosome==18,chr18.p$seg.gain.prop,NA)
      seg.gain.pro.chr18.p[length(seg.gain.pro.chr18.p)==0]=NA
      
      chr19.p=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="p",]
      seg.gain.pro.chr19.p=ifelse(chr19.p$chromosome==19,chr19.p$seg.gain.prop,NA)
      seg.gain.pro.chr19.p[length(seg.gain.pro.chr19.p)==0]=NA
      chr21.p=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="p",]
      seg.gain.pro.chr21.p=ifelse(chr21.p$chromosome==21,chr21.p$seg.gain.prop,NA)
      seg.gain.pro.chr21.p[length(seg.gain.pro.chr21.p)==0]=NA
      # . q arm
      chr3.q=hyper.s[hyper.s$chromosome==3 & hyper.s$arm=="q",]
      seg.gain.pro.chr3.q=ifelse(chr3.q$chromosome==3,chr3.q$seg.gain.prop,NA)
      seg.gain.pro.chr3.q[length(seg.gain.pro.chr3.q)==0]=NA
      chr5.q=hyper.s[hyper.s$chromosome==5 & hyper.s$arm=="q",]
      seg.gain.pro.chr5.q=ifelse(chr5.q$chromosome==5,chr5.q$seg.gain.prop,NA)
      seg.gain.pro.chr5.q[length(seg.gain.pro.chr5.q)==0]=NA
      chr7.q=hyper.s[hyper.s$chromosome==7 & hyper.s$arm=="q",]
      seg.gain.pro.chr7.q=ifelse(chr7.q$chromosome==7,chr7.q$seg.gain.prop,NA)
      seg.gain.pro.chr7.q[length(seg.gain.pro.chr7.q)==0]=NA
      chr9.q=hyper.s[hyper.s$chromosome==9 & hyper.s$arm=="q",]
      seg.gain.pro.chr9.q=ifelse(chr9.q$chromosome==9,chr9.q$seg.gain.prop,NA)
      seg.gain.pro.chr9.q[length(seg.gain.pro.chr9.q)==0]=NA
      chr11.q=hyper.s[hyper.s$chromosome==11 & hyper.s$arm=="q",]
      seg.gain.pro.chr11.q=ifelse(chr11.q$chromosome==11,chr11.q$seg.gain.prop,NA)
      seg.gain.pro.chr11.q[length(seg.gain.pro.chr11.q)==0]=NA
      chr15.q=hyper.s[hyper.s$chromosome==15 & hyper.s$arm=="q",]
      seg.gain.pro.chr15.q=ifelse(chr15.q$chromosome==15,chr15.q$seg.gain.prop,NA)
      seg.gain.pro.chr15.q[length(seg.gain.pro.chr15.q)==0]=NA
      
      chr18.q=hyper.s[hyper.s$chromosome==18 & hyper.s$arm=="q",]
      seg.gain.pro.chr18.q=ifelse(chr18.q$chromosome==18,chr18.q$seg.gain.prop,NA)
      seg.gain.pro.chr18.q[length(seg.gain.pro.chr18.q)==0]=NA
      
      chr19.q=hyper.s[hyper.s$chromosome==19 & hyper.s$arm=="q",]
      seg.gain.pro.chr19.q=ifelse(chr19.q$chromosome==19,chr19.q$seg.gain.prop,NA)
      seg.gain.pro.chr19.q[length(seg.gain.pro.chr19.q)==0]=NA
      chr21.q=hyper.s[hyper.s$chromosome==21 & hyper.s$arm=="q",]
      seg.gain.pro.chr21.q=ifelse(chr21.q$chromosome==21,chr21.q$seg.gain.prop,NA)
      seg.gain.pro.chr21.q[length(seg.gain.pro.chr21.q)==0]=NA
      
      perc.gain.chr3=unique(hyper.s[hyper.s$chromosome==3,"perc.gain.chr"])
      perc.gain.chr3[length(perc.gain.chr3)==0]=NA
      perc.gain.chr5=unique(hyper.s[hyper.s$chromosome==5,"perc.gain.chr"])
      perc.gain.chr5[length(perc.gain.chr5)==0]=NA
      perc.gain.chr7=unique(hyper.s[hyper.s$chromosome==7,"perc.gain.chr"])
      perc.gain.chr7[length(perc.gain.chr7)==0]=NA
      perc.gain.chr9=unique(hyper.s[hyper.s$chromosome==9,"perc.gain.chr"])
      perc.gain.chr9[length(perc.gain.chr9)==0]=NA
      perc.gain.chr11=unique(hyper.s[hyper.s$chromosome==11,"perc.gain.chr"])
      perc.gain.chr11[length(perc.gain.chr11)==0]=NA
      perc.gain.chr15=unique(hyper.s[hyper.s$chromosome==15,"perc.gain.chr"])
      perc.gain.chr15[length(perc.gain.chr15)==0]=NA
      
      perc.gain.chr18=unique(hyper.s[hyper.s$chromosome==18,"perc.gain.chr"])
      perc.gain.chr18[length(perc.gain.chr18)==0]=NA
      
      perc.gain.chr19=unique(hyper.s[hyper.s$chromosome==19,"perc.gain.chr"])
      perc.gain.chr19[length(perc.gain.chr19)==0]=NA
      perc.gain.chr21=unique(hyper.s[hyper.s$chromosome==21,"perc.gain.chr"])
      perc.gain.chr21[length(perc.gain.chr21)==0]=NA
      
      double.chrs=sum(hyper.s[hyper.s$chromosome %in% chroms.hyper,"arm.gain"])
      Cytogenetic="Trisomie"
      line.double.p=c(sam,type,whole.median.CN,seg.gain.pro.chr3.p,seg.gain.pro.chr5.p,seg.gain.pro.chr7.p,seg.gain.pro.chr9.p,
                      seg.gain.pro.chr11.p,seg.gain.pro.chr15.p,
                      seg.gain.pro.chr18.p,
                      seg.gain.pro.chr19.p,seg.gain.pro.chr21.p,double.chrs,Cytogenetic,"p",
                      perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,
                      perc.gain.chr18,
                      perc.gain.chr19,perc.gain.chr21)
      line.double.q=c(sam,type,whole.median.CN,seg.gain.pro.chr3.q,seg.gain.pro.chr5.q,seg.gain.pro.chr7.q,seg.gain.pro.chr9.q,
                      seg.gain.pro.chr11.q,seg.gain.pro.chr15.q,
                      seg.gain.pro.chr18.q,
                      seg.gain.pro.chr19.q,seg.gain.pro.chr21.q,double.chrs,Cytogenetic,"q",
                      perc.gain.chr3,perc.gain.chr5,perc.gain.chr7,perc.gain.chr9,perc.gain.chr11,perc.gain.chr15,
                      perc.gain.chr18,
                      perc.gain.chr19,perc.gain.chr21)
      hyper.1=rbind(hyper.1,line.double.p,line.double.q)
    }else{
      no.hyperdiploid.1=rbind(no.hyperdiploid.1,samples[iii])
    }
  }
  
  hyper1=as.data.frame(hyper.1)
  colnames(hyper1)=c("sample","type","whole.median.CN","seg.gain.pro.chr3","seg.gain.pro.chr5","seg.gain.pro.chr7","seg.gain.pro.chr9",
                     "seg.gain.pro.chr11","seg.gain.pro.chr15",
                     "seg.gain.pro.chr18",
                     "seg.gain.pro.chr19","seg.gain.pro.chr21","double.chrs","Cytogenetic","arm",
                     "perc.gain.chr3","perc.gain.chr5","perc.gain.chr7","perc.gain.chr9","perc.gain.chr11","perc.gain.chr15",
                     "perc.gain.chr18",
                     "perc.gain.chr19","perc.gain.chr21"
  )
  rownames(hyper1)=c(1:length(hyper1$sample))
  
  for(q in c(4:12,16:24)){
    hyper1[,q]=as.numeric(as.character(hyper1[,q]))
  }
  SAMPLE=unique(hyper1$sample)
  hyper2=NULL
  for(qq in 1:length(SAMPLE)){
    hyper1.sample=hyper1[hyper1$sample==SAMPLE[qq],]
    sam=SAMPLE[qq]
    type=unique(hyper1.sample$type)
    whole.median.CN=unique(hyper1.sample$whole.median.CN)
    cyto=unique(hyper1.sample$Cytogenetic)
    
    hyp.chr3=unique(ifelse(hyper1.sample$perc.gain.chr3>=0.6,1,0))
    hyp.chr5=unique(ifelse(hyper1.sample$perc.gain.chr5>=0.6,1,0))
    hyp.chr7=unique(ifelse(hyper1.sample$perc.gain.chr7>=0.6,1,0))
    hyp.chr9=unique(ifelse(hyper1.sample$perc.gain.chr9>=0.6,1,0))
    hyp.chr11=unique(ifelse(hyper1.sample$perc.gain.chr11>=0.6,1,0))
    hyp.chr15=unique(ifelse(hyper1.sample$perc.gain.chr15>=0.6,1,0))
    
    hyp.chr18=unique(ifelse(hyper1.sample$perc.gain.chr18>=0.6,1,0))
    
    hyp.chr19=unique(ifelse(hyper1.sample$perc.gain.chr19>=0.6,1,0))
    hyp.chr21=unique(ifelse(hyper1.sample$perc.gain.chr21>=0.6,1,0))
    
    hyp.sample=ifelse(sum(hyp.chr3,hyp.chr5,hyp.chr7,hyp.chr9,hyp.chr11,hyp.chr15,hyp.chr19,hyp.chr21,na.rm = T)>=2,"Hyperdiploidy","No.HRD") # Except chr18
    line=c(sam,type,whole.median.CN,
           hyp.chr3,hyp.chr5,hyp.chr7,hyp.chr9,hyp.chr11,hyp.chr15,hyp.chr19,hyp.chr21,hyp.chr18,
           hyp.sample,
           cyto)
    hyper2=rbind(line,hyper2)
  }
  hyper.2=as.data.frame(hyper2)
  colnames(hyper.2)=c("sample","type","whole.median.CN",
                      "hyp.chr3","hyp.chr5","hyp.chr7","hyp.chr9","hyp.chr11","hyp.chr15","hyp.chr19","hyp.chr21","hyp.chr18",
                      "hyp.sample","Cytogenetic")
  rownames(hyper.2)=c(1:length(hyper.2$sample))
  
  no.hyper=as.data.frame(no.hyperdiploid)
  if(length(no.hyper)>0){
    colnames(no.hyper)="sample"
  }
  no.hyper1=as.data.frame(no.hyperdiploid.1)
  if(length(no.hyper1)>0){
    colnames(no.hyper1)="sample"
  }
  
  no.hyper.1=rbind(no.hyper,no.hyper1)
  if(dim(no.hyper.1)[1]!=0){
    no.hyper.1$Cytogenetic="No_HRD"
  }
  
  hyper.5cohorts=hyper.2
  no.hyper.5cohorts=no.hyper.1
  
  hyper.def=hyper.5cohorts
  no.hyper.def=no.hyper.5cohorts
  
  write.table(
    hyper.def,
    paste(
      OUT.DIR,
      PERCENTAGE*100,
      "%_ARMS_",
      DATABASE,
      "hyperdiploid.txt",
      sep=""
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep="\t"
  )
  
  write.table(
    no.hyper.def,
    paste(
      OUT.DIR,
      PERCENTAGE*100,
      "%_ARMS_",
      DATABASE,
      "no.hyperdiploid.txt",
      sep=""
    ),
    col.names = T,
    row.names = F,
    quote = F,
    sep="\t"
  )
}




#### . . . - CytogeneticMatrixGenerator()

#### . . . - Function to create a cytogenetic matrix
####         FocalLargeMATRIX = 

CytogeneticMatrixGenerator=function(FocalLargeMATRIX,HY.FILE,DATABASE,OUT.DIR){
  
  ##################################################################################
  ##################################################################################
  ##############################################################################1###
  #                                                                                #
  #        .     .     .                                               .     .     .
  ####     .     .     .   PredictionModel - Cytogenetic abnormalies   .     .     .     ####
  #        .     .     .                                               .     .     .
  #                                                                                #
  ##################################################################################
  ##################################################################################
  ##################################################################################
  
  
  ###############################################################################################################################################################
  ######
  # . . . . .                                    PredictionModel - Matrix HRD chromosome
  ######
  ################################################################################################################################################################
  #####################
  ######## --- .                 
  #####################
  
  matrix=read.delim(FocalLargeMATRIX,
                    stringsAsFactors = F)
  pres.abs=matrix
  for(j in c(1:(ncol(matrix)-1))){
    pres.abs[,j]=ifelse(pres.abs[,j]>0,1,pres.abs[,j])
  }
  hyperdiploid=read.delim(HY.FILE,
                          stringsAsFactors = F)
  matrix1=merge(matrix,hyperdiploid,by.x="row.names",by.y="sample",all.x=T)
  pres.abs1=merge(pres.abs,hyperdiploid,by.x="row.names",by.y="sample",all.x=T)
  for(www in c(grep("hyp.chr",colnames(pres.abs1)))){
    pres.abs1[,www]=ifelse(is.na(pres.abs1[,www]),0,pres.abs1[,www])
  }
  pres.abs1$hyp.sample=ifelse(is.na(pres.abs1$hyp.sample),"No.HRD",pres.abs1$hyp.sample)
  matrix.hyper=pres.abs1[,c(1,grep("hyp.chr",colnames(pres.abs1)))]
  rownames(matrix.hyper)=matrix.hyper$Row.names
  matrix.hyper=matrix.hyper[,-1]
  colnames(matrix.hyper)=c("chr3+","chr5+","chr7+","chr9+","chr11+","chr15+","chr19+","chr21+","chr18+")
  write.table(matrix.hyper,
              paste(OUT.DIR,"/HRD_chromosomes",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
  #
  
  ##############################################################################################################################################################
  ######
  # . . . . .                                    PredictionModel - Matrix Other GISTIC peaks
  ######
  ################################################################################################################################################################
  #####################
  ######## --- .                 
  #####################
  
  matrix=read.delim(FocalLargeMATRIX,
                    stringsAsFactors = F)
  # gistic.amp=matrix[,c(4,5,8,16,17,18,20,24,25,27)]
  # for(j in c(1:(ncol(gistic.amp)))){
  #   gistic.amp[,j]=ifelse(gistic.amp[,j] %in% c(0.5,1),1,
  #                         ifelse(gistic.amp[,j] %in% c(1.5,2),2,
  #                                ifelse(gistic.amp[,j]==0,0,gistic.amp[,j])))
  # }
  gistic.amp=matrix[,c("Amp_14q11.2","Amp_14q32.33","Amp_17q22","Amp_2p11.2","Amp_2q24.3","Amp_2q32.1","Amp_4q13.2","Amp_6p24.3","Amp_6q11.1","Amp_8q24.21")]
  for(j in c(1:(ncol(gistic.amp)))){
    gistic.amp[,j]=ifelse(gistic.amp[,j] %in% c(0.5,1),1,
                          ifelse(gistic.amp[,j] %in% c(1.5,2),2,
                                 ifelse(gistic.amp[,j]==0,0,gistic.amp[,j])))
  }
  # gistic.del=matrix[,c(31:65)]
  # for(j in c(1:(ncol(gistic.del)-1))){
  #   gistic.del[,j]=ifelse(gistic.del[,j] %in% c(0.5,1),1,
  #                         ifelse(gistic.del[,j] %in% c(1.5,2),2,
  #                                ifelse(gistic.del[,j]==0,0,gistic.del[,j])))
  # }
  gistic.del=matrix[,c("Del_10p15.3","Del_10q24.32","Del_10q26.3","Del_11q22.1","Del_12p13.2","Del_12q24.31","Del_13q12.11","Del_13q14.2","Del_14q24.3","Del_14q32.33",
                       "Del_16p11.2","Del_16p13.3","Del_16q12.1","Del_16q22.2","Del_16q24.3","Del_17p13.1","Del_1p12","Del_1p22.1","Del_1p32.3","Del_20p13",   
                       "Del_20q13.12","Del_21p11.2","Del_22q11.22","Del_22q13.32","Del_2p11.2","Del_2q31.1","Del_2q37.3","Del_6q26","Del_7p11.2","Del_7p22.2",  
                       "Del_8p21.3","Del_8p23.1","Del_8q24.21","Del_9p21.3","database" )]
  for(j in c(1:(ncol(gistic.del)-1))){
    gistic.del[,j]=ifelse(gistic.del[,j] %in% c(0.5,1),1,
                          ifelse(gistic.del[,j] %in% c(1.5,2),2,
                                 ifelse(gistic.del[,j]==0,0,gistic.del[,j])))
  }
  
  amp=matrix[grep("Amp",colnames(matrix))]
  del=matrix[grep("Del",colnames(matrix))]
  same=intersect(gsub("Amp_","",colnames(amp)),gsub("Del_","",colnames(del)))
  matrix$sample=rownames(matrix)
  matrix[,c("sample",paste("Amp_",same,sep=""),paste("Del_",same,sep=""))]
  
  gistic.other=merge(gistic.amp,gistic.del,by="row.names")
  rownames(gistic.other)=gistic.other[,1]
  gistic.other=gistic.other[,-1]
  
  write.table(gistic.other,
              paste(OUT.DIR,"/gistic.other_All_peaks",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
  #
  
  
  ###############################################################################################################################################################
  ######
  # . . . . .                                    PredictionModel - Matrix 1q
  ######
  ################################################################################################################################################################
  #####################
  ######## --- .                 
  #####################
  
  matrix=read.delim(FocalLargeMATRIX,
                    stringsAsFactors = F)
  
  x1q=matrix[,c("Amp_1q21.1","Amp_1q21.3","Amp_1q44","database")]
  for(j in c(1:(ncol(x1q)-1))){
    x1q[,j]=ifelse(x1q[,j] %in% c(0.5,1),1,
                   ifelse(x1q[,j] %in% c(1.5,2),2,
                          ifelse(x1q[,j]==0,0,x1q[,j])))
  }
  
  write.table(x1q,
              paste(OUT.DIR,"/x1q_amp_gain",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
  #
  
  # ###############################################################################################################################################################
  # ######
  # # . . . . .                                    PredictionModel - Matrix Translocations
  # ######
  # ################################################################################################################################################################
  # #####################
  # ######## --- .                 
  # #####################
  # 
  # translocation=read.delim(paste(WF,"/input/translocations/commpass_translocations.txt",sep=""),
  #                          stringsAsFactors = F)
  # 
  # for(j in c(2:(ncol(translocation)))){
  #   translocation[,j]=ifelse(translocation[,j] %in% c(1),2,translocation[,j])
  # }
  # 
  # rownames(translocation)=translocation[,1]
  # translocation=translocation[,-1]
  # 
  # write.table(translocation,
  #             paste(WF,"/output/translocations.txt",sep=""),
  #             col.names = T,
  #             row.names = T,
  #             quote = F,
  #             sep = "\t")
  # #
}



#### . . . - MutationsMatrixGenerator()

#### . . . - Function to create the mutation matrix (Oncogenes, TSGs, Unknown Genes, X genes)
####         CODING.MUT = 
MutationsMatrixGenerator=function(CODING.MUT,METADATA,UCSC.gene,DRIVERs,GisticOtherAllPeaks,CNV,DATABASE,OUT.DIR){
  ##################################################################################
  ##################################################################################
  ##############################################################################1###
  #                                                                                #
  #        .     .     .                                               .     .     .
  ####     .     .     .   PredictionModel - Mutated Genes (Oncogenes, TSGs and Unknown)   .     .     .     ####
  #        .     .     .                                               .     .     .
  #                                                                                #
  ##################################################################################
  ##################################################################################
  ##################################################################################
  ##############################################################################################################################################################
  ######
  # . . . . .                                    PredictionModel - ONCOGENE
  ######
  ################################################################################################################################################################
  
  #####################
  ######## --- .       Check old matrix          
  #####################
  
  # matrix=read.delim(paste(WF,"/input/commpass_info_test.txt",sep=""),
  #                   stringsAsFactors = F)
  
  UCSC.HG19=UCSC.gene
  
  #####################
  ########
  ######## --- .       Selection drivers            ######
  ########
  #####################
  drivers=read.delim(DRIVERs,
                     stringsAsFactors = F)
  drivers$n.CallBy.new=drivers$DNDS+drivers$MUTSIG+drivers$ONCODRIVER+drivers$FISHHOOK
  drivers1=drivers[drivers$n.CallBy.new>1 | drivers$MM.driver=="MM",]
  colnames(drivers1)[1]="Gene.Symbol"
  unique(drivers1[order(drivers1$Gene.Symbol),"Gene.Symbol"])
  drivers2=drivers1[,c("Gene.Symbol","Role","Role_in_Cancer","MM.driver")]
  
  #####################
  ########
  ######## --- .       Oncogenes             ######        
  ########
  #####################
  oncogene=drivers2[drivers2$Role=="Oncogene",]
  oncogene=oncogene[complete.cases(oncogene$Gene.Symbol),]
  
  #########################################################################################################
  ########
  ######## --- .                                    Cohort                      ######
  ########
  #########################################################################################################
  # . . . - All mutation file by DNDS analysis
  mut.dnds=read.delim(CODING.MUT,
                      stringsAsFactors = F)
  mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
  mut.dnds.gene=mut.dnds[mut.dnds$gene %in% oncogene$Gene.Symbol,]
  
  mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
  mat=as.matrix(mat)
  mat=ifelse(mat>1,1,mat)
  mat=as.data.frame(mat)
  
  matrix=read.delim(METADATA,
                    stringsAsFactors = F)
  samples.NA=matrix[matrix$snv==0,"sample"]
  #- Using the old matrix - NA matrix
  matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
  colnames(matrix.NA)=colnames(mat)
  rownames(matrix.NA)=samples.NA
  # Matrix with the samples without mutations
  samples.0=setdiff(matrix[matrix$snv==1,"sample"],c(rownames(mat),rownames(matrix.NA)))
  matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
  colnames(matrix.0)=colnames(mat)
  rownames(matrix.0)=samples.0
  
  matrix.oncogene=rbind(mat,matrix.0,matrix.NA)

  dim(matrix.oncogene)
  
  missing=setdiff(unique(oncogene$Gene.Symbol),colnames(matrix.oncogene))
  for(i in 1:length(missing)){
    matrix.oncogene[,ncol(matrix.oncogene)+1]=0
    matrix.oncogene[,ncol(matrix.oncogene)]=ifelse(rownames(matrix.oncogene) %in% samples.NA,NA,matrix.oncogene[,ncol(matrix.oncogene)])
    colnames(matrix.oncogene)[ncol(matrix.oncogene)]=missing[i]
  }
  matrix.oncogene=matrix.oncogene[,order(colnames(matrix.oncogene))]
  
  write.table(matrix.oncogene,
              paste(OUT.DIR,"/oncogene_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep="\t")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #
  # . . . -> REMEMBER: - - - For ATRX,KDM5C,KDM6A,RPL10 we consider only the mutations because they are in the X chromosome.
  #
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  #####################
  ########
  ######## --- .       Tumor Suppressor Genes - TSG         . ---########        
  ########
  #####################
  tsg=drivers2[drivers2$Role=="TSG",]
  tsg=tsg[complete.cases(tsg$Gene.Symbol),]
  
  #########################################################################################################
  ########
  ######## --- .                                    Cohort                      ######
  ########
  #########################################################################################################
  # . . . - All mutation file by DNDS analysis
  mut.dnds=read.delim(CODING.MUT,
                      stringsAsFactors = F)
  mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
  mut.dnds.gene=mut.dnds[mut.dnds$gene %in% tsg$Gene.Symbol,]
  
  mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
  mat=as.matrix(mat)
  mat=ifelse(mat>1,1,mat)
  mat=as.data.frame(mat)
  samples.NA=matrix[matrix$snv==0,"sample"]
  # -  NA matrix
  matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
  colnames(matrix.NA)=colnames(mat)
  rownames(matrix.NA)=samples.NA
  # - Matrix with the samples without mutations
  samples.0=setdiff(matrix[matrix$snv==1,"sample"],c(rownames(mat),rownames(matrix.NA)))
  matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
  colnames(matrix.0)=colnames(mat)
  rownames(matrix.0)=samples.0
  
  matrix.tsg=rbind(mat,matrix.0,matrix.NA)
  
  missing=setdiff(unique(tsg$Gene.Symbol),colnames(matrix.tsg))
  for(i in 1:length(missing)){
    matrix.tsg[,ncol(matrix.tsg)+1]=0
    matrix.tsg[,ncol(matrix.tsg)]=ifelse(rownames(matrix.tsg) %in% samples.NA,NA,matrix.tsg[,ncol(matrix.tsg)])
    colnames(matrix.tsg)[ncol(matrix.tsg)]=missing[i]
  }
  matrix.tsg=matrix.tsg[,order(colnames(matrix.tsg))]
  
  write.table(matrix.tsg,
              paste(OUT.DIR,"/tsg_mutated_genes_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep="\t")
  
  
  #####
  # . . . - Upload GISTIC file to performe a comparison at the end of the analysis
  gistic=read.delim(GisticOtherAllPeaks,
                    stringsAsFactors = F)
  gistic=gistic[,grep("Del_",colnames(gistic))]
  # rownames(gistic)=gsub("_1_BM","",rownames(gistic))
  # gistic.study=gistic[rownames(gistic) %in% matrix[matrix$study=="MMRF","sample"],]
  
  gistic.study.1=gistic[,c("Del_1p12","Del_17p13.1","Del_1p32.3",
                           "Del_1p22.1"),]
  colnames(gistic.study.1)=c("Del_TENT5C","Del_TP53","Del_CDKN2C",
                             "Del_RPL5")
  #####
  
  #######
  # . . . - Upload study CNV file to check the gene region CNV values 
  #######
  cnv=CNV
  # cnv$sample=gsub("_1_BM","",cnv$sample)
  
  tsg$Gene.Symbol
  
  ref.genes=read.delim(UCSC.HG19,
                       stringsAsFactors = F)
  ts.genes=ref.genes[ref.genes$geneSymbol %in% tsg$Gene.Symbol,]
  ts.genes=ts.genes[,c("chrom","txStart","txEnd",
                       "geneSymbol")]
  ts.genes=ts.genes %>%
    distinct(geneSymbol,.keep_all = T)
  cnv.1=cnv
  # length(unique(cnv.1$sample))
  # head(cnv.1)
  # head(ts.genes)
  colnames(ts.genes)=c("chrom","start","end","Symbol")
  ts.genes$chrom=gsub("chr","",ts.genes$chrom)
  
  
  gr.cnv=with(cnv.1, GRanges(chrom, IRanges(start=start-50000, end=end+50000)))
  values(gr.cnv)=cnv.1[,c(1,5:ncol(cnv.1))]
  
  gr.tsg=with(ts.genes,GRanges(chrom,IRanges(start=start,end=end)))
  values(gr.tsg)=ts.genes[,c(4:ncol(ts.genes))]
  
  ranges <- merge(as.data.frame(gr.cnv),
                  as.data.frame(gr.tsg),
                  by="seqnames",
                  suffixes=c("A","B"))
  
  ranges_1=ranges[with(ranges, startB <= startA & endB <= endA & endB > startA),]
  ranges_2=ranges[with(ranges, startB >= startA & endB >= endA & startB < endA),]
  ranges_3=ranges[with(ranges, startB >= startA & endB <= endA),]
  ranges_4=ranges[with(ranges, startB < startA & endB > endA),]
  ranges.all<- rbind.data.frame(ranges_1, ranges_2, ranges_3, ranges_4)
  
  ranges.all=unique(ranges.all)
  ranges.all.1=ranges.all %>%
    arrange(sample,seqnames,startA)
  ranges.all.2=ranges.all.1[,c("sample","seqnames","startA","endA","total","minor","X")]
  # unique(ranges.all.2$total)
  # table(ranges.all.2$total)
  # table(round(ranges.all.2$total,0))
  
  ################################################################################
  # TOTAL CNV -> <1.5 will be loss, 2.5 will be consider diploid, 3.5 a triploid and 4.5 tetraploid
  # MINOR CNV -> <0.5 will be homozygous deletion (0)
  
  ranges.all.2$total=ifelse(ranges.all.2$total<2.5 & ranges.all.2$total>=1.5,2,
                            ifelse(ranges.all.2$total<1.5 & ranges.all.2$total>=0.5,1,
                                   ifelse(ranges.all.2$total<3.5 & ranges.all.2$total>=2.5,3,
                                          ifelse(ranges.all.2$total<4.5 & ranges.all.2$total>=3.5,4,ranges.all.2$total))))
  
  
  
  ranges.all.2$minor=ifelse(ranges.all.2$minor<0.5,0,ranges.all.2$minor)
  ################################################################################
  
  
  ranges.all.2$major=ranges.all.2$total-ranges.all.2$minor
  ranges.all.2$code=ifelse(ranges.all.2$major+ranges.all.2$minor>=1.5 & ranges.all.2$minor>=0.5,0,
                           ifelse(ranges.all.2$major+ranges.all.2$minor<1.5 | (ranges.all.2$minor<0.5 & ranges.all.2$major>=1.5),1,
                                  ifelse(ranges.all.2$total<0.5,2,NA)))
  
  
  ranges.all.2$code=ifelse(is.na(ranges.all.2$code) & ranges.all.2$total<1.5,1,
                           ifelse(is.na(ranges.all.2$code) & ranges.all.2$total<0.5,2,
                                  ifelse(is.na(ranges.all.2$code) & ranges.all.2$total>=1.5 & ranges.all.2$total<2.5,0,
                                         ifelse(is.na(ranges.all.2$code) & ranges.all.2$total>=2.5,0,ranges.all.2$code)
                                  )
                           )
  )
  
  ranges.all.2=setDT(ranges.all.2)
  matrix.cnv=data.table::dcast(unique(ranges.all.2[,c("X","sample","code")]),
                               formula=sample~X,
                               toString,
                               value.var = "code")
  #data.table::setDF(matrix.cnv)
  matrix.cnv=as.data.frame(matrix.cnv)
  rownames(matrix.cnv)=matrix.cnv$sample
  matrix.cnv=matrix.cnv[,-1]
  for(j in 1:ncol(matrix.cnv)){
    matrix.cnv[,j]=ifelse(matrix.cnv[j]=="",0,matrix.cnv[,j])
  }
  for(i in 1:length(colnames(matrix.cnv))){
    for(z in 1:length(rownames(matrix.cnv))){
      cel=matrix.cnv[z,i]
      c1=str_split_fixed(cel,",",Inf)
      c1=as.matrix(c1)
      c1=gsub(" ","",c1)
      c1=as.data.frame(c1)
      for(zz in 1:length(colnames(c1))){
        c1[,zz]=as.numeric(c1[,zz])
      }
      c1$max=max(c1)
      matrix.cnv[z,i]=c1$max
    }
  }
  
  
  samples.NA=matrix[matrix$cnv==0,"sample"]
  # -  NA matrix
  matrix.NA=matrix(NA,ncol=ncol(matrix.cnv),nrow=length(samples.NA))
  colnames(matrix.NA)=colnames(matrix.cnv)
  rownames(matrix.NA)=samples.NA
  # - Matrix with the samples without mutations
  samples.0=setdiff(cnv$sample,c(rownames(matrix.cnv),rownames(matrix.NA)))
  matrix.0=matrix(0,ncol=ncol(matrix.cnv),nrow=length(samples.0))
  colnames(matrix.0)=colnames(matrix.cnv)
  rownames(matrix.0)=samples.0
  
  matrix.tsg=rbind(matrix.cnv,matrix.0,matrix.NA)
  dim(matrix.tsg)
  
  colnames(matrix.tsg)=paste("Del_",colnames(matrix.tsg),sep="")
  
  
  write.table(matrix.tsg,
              paste(OUT.DIR,"/tsg_deleted_genes_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep="\t")
  
  ##################################################################################################################################
  ##################################################################################################################################
  
  # . . . Merge two matrices
  # Mutated gene matrix
  mut=read.delim(paste(OUT.DIR,"/tsg_mutated_genes_",DATABASE,".txt",sep=""),
                 stringsAsFactors = F)
  missing1=setdiff(tsg$Gene.Symbol,colnames(mut))
  if(length(missing1)>=1){
    for(i in 1:length(missing1)){
      mut[,ncol(mut)+1]=0
      colnames(mut)[ncol(mut)]=missing1[i]
    }
  }
  
  mut=mut[,order(colnames(mut))]
  
  colnames(mut)=paste("Mut_",colnames(mut),sep="")
  
  
  # Deleted gene matrix
  cnv=read.delim(paste(OUT.DIR,"/tsg_deleted_genes_",DATABASE,".txt",sep=""),
                 stringsAsFactors = F)
  
  missing2=setdiff(tsg$Gene.Symbol,gsub("Del_","",colnames(cnv)))
  if(length(missing2)>=1){
    for(i in 1:length(missing2)){
      cnv[,ncol(cnv)+1]=0
      colnames(cnv)[ncol(cnv)]=paste("Del_",missing2[i],sep="")
    }
  }
  
  cnv=cnv[,order(colnames(cnv))]
  
  study.matrix=merge(cnv,mut,by="row.names",all=T)
  
  study.matrix.1=study.matrix
  rownames(study.matrix.1)=study.matrix.1$Row.names
  study.matrix.1=study.matrix.1[,-1]
  
  for(i in 1:ncol(study.matrix.1)){
    study.matrix.1[,i]=ifelse(is.na(study.matrix.1[,i]),999,study.matrix.1[,i])
  }
  
  for(j in 1:36){
    study.matrix.1[,72+j]=rowSums(study.matrix.1[,c(j,36+j)],na.rm = T)
    colnames(study.matrix.1)[ncol(study.matrix.1)]=gsub("Del_","",colnames(study.matrix.1)[j])
  }
  
  for(z in 1:72){
    study.matrix.1[,z]=ifelse(study.matrix.1[,z]==999,NA,study.matrix.1[,z])
  }
  for(z in 73:ncol(study.matrix.1)){
    study.matrix.1[,z]=ifelse(study.matrix.1[,z] %in% c(3,1001),2,
                              ifelse(study.matrix.1[,z]==1000,1,
                                     ifelse(study.matrix.1[,z]==999,0,
                                            ifelse(study.matrix.1[,z]==1998,NA,study.matrix.1[,z]))))
  }
  
  # ATRX,KDM5C,KDM6A,RPL10
  study.matrix.2=study.matrix.1[,c(73:108)]
  colnames(study.matrix.2)=gsub("Mut_","",colnames(study.matrix.2))
  study.matrix.2=study.matrix.2[,order(colnames(study.matrix.2))]
  
  write.table(study.matrix.2,
              paste(OUT.DIR,"/tsg_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep="\t")
  
  
  
  #####################
  ########
  ######## --- .       Unknown         . ---########        
  ########
  #####################
  ukn=drivers2[drivers2$Role=="Unknown",]
  ukn1=drivers2[is.na(drivers2$Role),]
  
  ukn=rbind(ukn,ukn1)
  ukn=ukn[complete.cases(ukn$Gene.Symbol),]
  
  
  #########################################################################################################
  ########
  ######## --- .                                    Cohort                                   ######
  ########
  #########################################################################################################
  # . . . - All mutation file by DNDS analysis
  mut.dnds=read.delim(CODING.MUT,
                      stringsAsFactors = F)
  mut.dnds=mut.dnds[mut.dnds$impact != "Synonymous",]
  mut.dnds.gene=mut.dnds[mut.dnds$gene %in% ukn$Gene.Symbol,]
  
  mat=as.data.frame.matrix(table(mut.dnds.gene$sampleID,mut.dnds.gene$gene))
  mat=as.matrix(mat)
  mat=ifelse(mat>1,1,mat)
  mat=as.data.frame(mat)
  samples.NA=matrix[matrix$snv==0,"sample"]
  # -  NA matrix
  matrix.NA=matrix(NA,ncol=ncol(mat),nrow=length(samples.NA))
  colnames(matrix.NA)=colnames(mat)
  rownames(matrix.NA)=samples.NA
  # - Matrix with the samples without mutations
  samples.0=setdiff(matrix[matrix$snv==1,"sample"],c(rownames(mat),rownames(matrix.NA)))
  matrix.0=matrix(0,ncol=ncol(mat),nrow=length(samples.0))
  colnames(matrix.0)=colnames(mat)
  rownames(matrix.0)=samples.0
  
  matrix.ukn=rbind(mat,matrix.0,matrix.NA)
  dim(matrix.ukn)
  
  missing=setdiff(unique(ukn$Gene.Symbol),colnames(matrix.ukn))
  for(i in 1:length(missing)){
    matrix.ukn[,ncol(matrix.ukn)+1]=0
    matrix.ukn[,ncol(matrix.ukn)]=ifelse(rownames(matrix.ukn) %in% samples.NA,NA,matrix.ukn[,ncol(matrix.ukn)])
    colnames(matrix.ukn)[ncol(matrix.ukn)]=missing[i]
  }
  matrix.ukn=matrix.ukn[,order(colnames(matrix.ukn))]
  
  
  write.table(matrix.ukn,
              paste(OUT.DIR,"/unk_mutated_genes_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep="\t")
  
  
  #########################################################################################################
  ########
  ######## --- .                            Chromosome X Genes                               . --- ########   
  ########
  #########################################################################################################
  
  # "ATRX","KDM5C","KDM6A","RPL10"
  
  commpass=read.delim(paste(OUT.DIR,"/tsg_mutated_genes_",DATABASE,".txt",sep=""),
                      stringsAsFactors = F)
  commpass=commpass[,colnames(commpass) %in% c("ATRX","KDM5C","KDM6A","RPL10")]
  
  write.table(commpass,
              paste(OUT.DIR,"/chrx_tsg_",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
}







#### . . . - APOBECanalyzer()

#### . . . - Function to create the mutation matrix (Oncogenes, TSGs, Unknown Genes, X genes)
####         SNV.OBJECT = 
APOBECanalyzer=function(SNV.OBJECT,ITER=10,SAMPLE,CHROM,POS,REF,ALT,DATABASE,OUT.DIR){
  
  source("~/Google Drive/My Drive/Bachisio/Desktop_OCT2023/Tools/R_functions/mmsig-master/R/mmsig_functions_hg38.R")
  ##################################################################################
  ##################################################################################
  ##############################################################################1###
  #                                                                                #
  #        .     .     .                                               .     .     .
  ####     .     .     .   PredictionModel - Mutational Signatures   .     .     .     ####
  #        .     .     .                                               .     .     .
  #                                                                                #
  ##################################################################################
  ##################################################################################
  ##################################################################################
  set.seed(1)
  
  rotatedAxisElementText = function(angle,position='x'){
    angle     = angle[1];
    position  = position[1]
    positions = list(x=0,y=90,top=180,right=270)
    if(!position %in% names(positions))
      stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
    if(!is.numeric(angle))
      stop("'angle' must be numeric",call.=FALSE)
    rads  = (angle - positions[[ position ]])*pi/180
    hjust = 0.5*(1 - sin(rads))
    vjust = 0.5*(1 + cos(rads))
    element_text(angle=angle,vjust=vjust,hjust=hjust)
  }
  
  scale_fill_sigs <- function(...){
    ggplot2:::manual_scale('fill',
                           values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"),"lightpink",
                                               "firebrick","lightblue1","darkslategray4"),
                                             c("SBS1", "SBS2", "SBS5","SBS8", "SBS9","SBS13","SBS18","SBS-MM1",
                                               "SBS17b","SBS31","SBS35","E_SBS37")),
                           ...)
  }
  
  # . . . - Signatures definition
  sig_ref <- signature_ref[c("class","sub", "tri", "SBS1", "SBS2", "SBS5", "SBS8", 
                             "SBS9", "SBS13", "SBS18")]
  sig<-sig_ref[,-1]
  sig<-sig[,c("sub", "tri","SBS1","SBS2","SBS13","SBS5", "SBS8", 
              "SBS9", "SBS18")]
  ##################################################################################
  ##################################################################################

  ########################################################################################################################
  
  dnds.input=SNV.OBJECT
  
  dnds.input=dnds.input[,c(SAMPLE,CHROM,POS,REF,ALT)]
  colnames(dnds.input)=c("sample","chr","pos","ref","alt")
  
  dnds.input=dnds.input[(nchar(dnds.input$ref)+nchar(dnds.input$alt))==2,]
  dnds.input$chr=gsub("chr","",dnds.input$chr)
  dnds.input$chr=paste("chr",dnds.input$chr,sep="")
  
  
  sig_out100 <- mm_fit_signatures(muts.input=dnds.input,
                                  sig.input=sig,
                                  input.format = "vcf",
                                  strandbias = T,
                                  bootstrap = TRUE,
                                  iterations = ITER,
                                  refcheck=FALSE,
                                  cos_sim_threshold = 0.01,
                                  force_include = c("SBS1","SBS5"),
                                  dbg=FALSE)
  sig_order = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9","SBS17b", 
                "SBS18","SBS-MM1","SBS31","E_SBS25","SBS35","E_SBS37")
  est<- sig_out100$estimate
  
  pr=est[est$SBS2!=0 | est$SBS13!=0,]
  bootSigsPlot(sig_out100$bootstrap)
  
  apobec=sig_out100$bootstrap[sig_out100$bootstrap$signature %in% c("SBS2","SBS13"),]
  
  samples=unique(apobec$sample)
  APOBEC=NULL
  for(i in 1:length(samples)){
    sub=apobec[apobec$sample==samples[i],]
    estimate=sub[1,6]+sub[2,6]
    row=c(samples[i],estimate)
    APOBEC=rbind(APOBEC,row)
  }
  APOBEC=as.data.frame(APOBEC)
  rownames(APOBEC)=1:length(APOBEC$V1)
  colnames(APOBEC)=c("sample","apobec")
  APOBEC$apobec=as.numeric(as.character(APOBEC$apobec))
  
  summary(APOBEC$apobec)
  
  
  APOBEC$offset=0
  rownames(APOBEC)=APOBEC$sample
  APOBEC=APOBEC[,-1]
  
  write.table(APOBEC,
              paste(OUT.DIR,"APOBEC_relative_contribution_matrix",DATABASE,".txt",sep=""),
              col.names = T,
              row.names = T,
              quote = F,
              sep = "\t")
}






#### . . . - MatrixGenerator()

#### . . . - Function to create the mutation matrix (Oncogenes, TSGs, Unknown Genes, X genes)
####         INPUT.FILE.FOLDER = 
FinalMatrixGenerator=function(INPUT.FILE.FOLDER,DATABASE,OUTPUT){
  ##################################################################################
  ##################################################################################
  ##############################################################################1###
  #                                                                                #
  #        .     .     .                                               .     .     .
  ####     .     .     .   PredictionModel - Matrix Generation   .     .     .     ####
  #        .     .     .                                               .     .     .
  #                                                                                #
  ##################################################################################
  ##################################################################################
  
  ##################################################################################
  ##################################################################################
  ##################################################################################
  ##################################################################################
  
  
  
  ##############################################################################################################################################################
  ######
  # . . . . .                                    PredictionModel - Matrix Creation 
  ######
  ################################################################################################################################################################
  
  #######
  #### . . .    - HRD Chromosomes -     . . . ####
  #######
  hrd=read.delim(paste(INPUT.FILE.FOLDER,"HRD_chromosomes",DATABASE,".txt",sep=""),stringsAsFactors = F)
  colnames(hrd)=paste(colnames(hrd),"gain",sep="")
  colnames(hrd)=paste("CNV_",colnames(hrd),sep="")
  
  #######
  #### . . .    - 1q -     . . . ####
  #######
  x1q=read.delim(paste(INPUT.FILE.FOLDER,"x1q_amp_gain",DATABASE,".txt",sep=""),stringsAsFactors = F)
  colnames(x1q)=paste("CNV_",colnames(x1q),sep="")
  # . - Merge CLINIC info with 1q info 
  clinic1=merge(hrd,x1q,by="row.names",all=T)
  clinic1=clinic1[,-ncol(clinic1)]
  colnames(clinic1)[1]="sample"
  
  #######
  #### . . .    - GISTIC peaks -     . . . ####
  #######
  gistic=read.delim(paste(INPUT.FILE.FOLDER,"gistic.other_All_peaks",DATABASE,".txt",sep=""),stringsAsFactors = F)
  # Remove the GISTIC peaks with TSGs
  # Del_1p32.3,Del_1p12,Del_1p22.1,Del_13q14.2,Del_16p13.3,Del_16q12.1,Del_17p13.1
  gistic=gistic[,!colnames(gistic) %in% c("Del_1p32.3","Del_1p12","Del_1p22.1","Del_13q14.2","Del_16p13.3","Del_16q12.1","Del_17p13.1",
                                          "Amp_14q32.33","Del_14q32.33","Del_22q11.22","Del_2p11.2","Amp_2p11.2")]
  gistic=gistic[,-ncol(gistic)]
  for(i in 1:ncol(gistic)){
    gistic[,i]=ifelse(gistic[,i]>1,1,gistic[,i])
  }
  colnames(gistic)=paste("CNV_",colnames(gistic),sep="")
  # . - Merge CLINIC info with HRD info 
  clinic2=merge(clinic1,gistic,by.x="sample","by.y"="row.names",all=T)
  
  #######
  #### . . .    - ONCOGENEs -     . . . ####
  #######
  oncogenes=read.delim(paste(INPUT.FILE.FOLDER,"oncogene_",DATABASE,".txt",sep=""),stringsAsFactors = F)
  for(i in 1:ncol(oncogenes)){
    oncogenes[,i]=ifelse(oncogenes[,i]>1,1,oncogenes[,i])
  }
  colnames(oncogenes)=paste("SNV_",colnames(oncogenes),sep="")
  # . - Merge CLINIC info with HRD info 
  clinic3=merge(clinic2,oncogenes,by.x="sample","by.y"="row.names",all=T)
  
  #######
  #### . . .    - TSGs -     . . . ####
  #######
  tsg=read.delim(paste(INPUT.FILE.FOLDER,"tsg_",DATABASE,".txt",sep=""),stringsAsFactors = F)
  colnames(tsg)=paste("CNV.SNV_",colnames(tsg),sep="")
  # . - Merge CLINIC info with HRD info 
  clinic4=merge(clinic3,tsg,by.x="sample","by.y"="row.names",all=T)
  
  #######
  #### . . .    - UNKNOWNs -     . . . ####
  #######
  unk=read.delim(paste(INPUT.FILE.FOLDER,"unk_mutated_genes_",DATABASE,".txt",sep=""),stringsAsFactors = F)
  colnames(unk)=paste("SNV_",colnames(unk),sep="")
  # . - Merge CLINIC info with HRD info 
  clinic5=merge(clinic4,unk,by.x="sample","by.y"="row.names",all=T)
  
  #######
  #### . . .    - APOBEC -     . . . ####
  #######
  apobec=read.delim(paste(INPUT.FILE.FOLDER,"APOBEC_relative_contribution_matrix",DATABASE,".txt",sep=""),stringsAsFactors = F)
  # . - Merge CLINIC info with APOBEC info 
  clinic6=merge(clinic5,apobec,by.x="sample","by.y"="row.names",all=T)
  # remove offset column
  clinic6=clinic6[,-ncol(clinic6)]
  
  #######
  #### . . .    - APOBEC presence or absence -     . . . ####
  #######
  #deciles=quantile(na.omit(clinic6$apobec), probs = seq(0, 1, 1/10),)
  # . . . - These values are the deciles obtained using the whole cohort.
  #         0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
  # 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.02891755 0.04139446 0.06430805 0.11203844 0.91384142 
  # clinic6$APOBEC=ifelse(clinic6$apobec<deciles[7],0,
  #                       ifelse(clinic6$apobec>=deciles[7] & clinic6$apobec<deciles[10],1,
  #                              ifelse(clinic6$apobec>=deciles[10],2,NA)))
  clinic6$APOBEC=ifelse(clinic6$apobec<0.02891755,0,
                        ifelse(clinic6$apobec>=0.02891755 & clinic6$apobec<0.11203844,1,
                               ifelse(clinic6$apobec>=0.11203844,2,NA)))
  
  #######
  #### . . .    - ChrX genes -     . . . ####
  #######
  clinic6=clinic6[,!colnames(clinic6) %in% c("CNV.SNV_ATRX","CNV.SNV_KDM5C","CNV.SNV_KDM6A","CNV.SNV_RPL10")]
  chrX=read.delim(paste(INPUT.FILE.FOLDER,"chrx_tsg_",DATABASE,".txt",sep=""),
                  stringsAsFactors = F)
  colnames(chrX)=paste("SNV_",colnames(chrX),sep="")
  clinic7=merge(clinic6,chrX,by.x="sample",by.y="row.names",all.x=T)
  
  #######
  #### . . .    - Matrix -     . . . ####
  #######
  
  clinic7=clinic7[,-grep("IGLL5",colnames(clinic7))] # IGLL5
  
  clinic8=clinic7[,c(1, # Clinic
                     2:45, # copy number
                     68:99, # TSGs
                     46:67,100:132,135:138, # Oncogenes and Unknown and chrX
                     134 # APOBEC
  )]
  
  
  clinic9=clinic8
  #### Collapse 1q peaks - . . . ####
  clinic9$Gain_Amp1q=paste(clinic9$CNV_Amp_1q21.1,clinic9$CNV_Amp_1q21.3,clinic9$CNV_Amp_1q44,sep="_")
  #table(clinic9$Gain_Amp1q)
  clinic9$Gain_Amp1q=ifelse(clinic9$Gain_Amp1q=="0_0_0",0,
                            ifelse(clinic9$Gain_Amp1q=="NA_NA_NA",NA,
                                   ifelse(clinic9$Gain_Amp1q %in% c("1_1_1","1_1_0","0_1_0","0_1_1","0_0_1","1_0_0","1_0_1",
                                                                    "NA_1_0","NA_1_1","NA_NA_1","1_NA_NA","NA_1_NA"),1,
                                          ifelse(clinic9$Gain_Amp1q %in% c("1_2_1","0_2_0","2_1_1","0_2_1","1_2_0","2_2_0","2_2_2","2_2_1","0_1_2",
                                                                           "2_1_0","1_1_2","1_2_2","0_0_2","0_2_2","2_0_0",
                                                                           "2_NA_NA","NA_2_2","NA_1_2","NA_2_1","NA_2_NA","NA_NA_2","2_NA_1","2_NA_2","1_NA_2",
                                                                           "1_2_NA","2_1_NA","2_2_NA"),2,NA))))
  colnames(clinic9)[ncol(clinic9)]="CNV_Gain_Amp1q"
  #### . . . - Collapse 8p - . . . ####
  clinic9$CNV_chr8p.loss=paste(clinic9$CNV_Del_8p21.3,clinic9$CNV_Del_8p23.1,sep="_")
  #table(clinic9$CNV_chr8p.loss,useNA="ifany")
  clinic9$CNV_chr8p.loss=ifelse(clinic9$CNV_chr8p.loss=="0_0",0,
                                ifelse(clinic9$CNV_chr8p.loss=="NA_NA",NA,
                                       ifelse(clinic9$CNV_chr8p.loss %in% c("1_1","1_0","0_1"),1,NA)))
  
  #### . . . - Column selection - . . . ####
  # clinic10=clinic9[,c("sample",
  #                     "CNV_chr3.gain","CNV_chr5.gain",   
  #                     "CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_chr18.gain",  
  #                     "CNV_Gain_Amp1q","CNV_Amp_14q11.2","CNV_Amp_17q22","CNV_Amp_2q24.3","CNV_Amp_2q32.1","CNV_Amp_4q13.2","CNV_Amp_6p24.3",  
  #                     "CNV_Amp_8q24.21","CNV_chr8p.loss","CNV_Del_10p15.3","CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_11q22.1","CNV_Del_12p13.2", 
  #                     "CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_20p13","CNV_Del_20q13.12","CNV_Del_21p11.2","CNV_Del_22q13.32","CNV_Del_2q31.1",  
  #                     "CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_7p11.2","CNV_Del_7p22.2","CNV_Del_8q24.21","CNV_Del_9p21.3","CNV.SNV_ARID1A",  
  #                     "CNV.SNV_ARID2","CNV.SNV_ATM","CNV.SNV_BTG1","CNV.SNV_CDKN1B","CNV.SNV_CDKN2C","CNV.SNV_CREBBP","CNV.SNV_CYLD",    
  #                     "CNV.SNV_DNMT3A","CNV.SNV_EP300","CNV.SNV_FAM46C","CNV.SNV_FUBP1","CNV.SNV_KMT2B","CNV.SNV_KMT2C","CNV.SNV_MAX",     
  #                     "CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_NFKBIA","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_RASA2","CNV.SNV_RB1",     
  #                     "CNV.SNV_RPL5","CNV.SNV_SETD2","CNV.SNV_SP140","CNV.SNV_TET2","CNV.SNV_TGDS","CNV.SNV_TP53","CNV.SNV_TRAF2",   
  #                     "CNV.SNV_TRAF3","CNV.SNV_UBR5","CNV.SNV_ZNF292","SNV_BRAF","SNV_CCND1","SNV_DIS3","SNV_FGFR3",       
  #                     "SNV_H1.4","SNV_IDH1","SNV_IDH2","SNV_IRF4","SNV_KRAS","SNV_NFKB2","SNV_NRAS",        
  #                     "SNV_PIK3CA","SNV_PTPN11","SNV_SF3B1","SNV_ABCF1","SNV_ACTG1","SNV_BHLHE41","SNV_DDX3X",       
  #                     "SNV_DTX1","SNV_DUSP2","SNV_H1.5","SNV_H1.2","SNV_H1.3","SNV_H2BC12","SNV_HLA.C",       
  #                     "SNV_HUWE1","SNV_INO80E","SNV_IRF1","SNV_KLHL6","SNV_LCE1D","SNV_LTB","SNV_MAN2C1",      
  #                     "SNV_MRGPRX4","SNV_NDUFC2","SNV_PABPC1","SNV_PANK3","SNV_PRKD2","SNV_RFTN1","SNV_RPRD1B",      
  #                     "SNV_RPS3A","SNV_SAMHD1","SNV_SLC35G5","SNV_SYT17","SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7",        
  #                     "SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_BCL7A","SNV_EGR1",        
  #                     "SNV_MAF","SNV_MAFB","SNV_MAML2","SNV_PIM1","SNV_TCL1A","SNV_XBP1",         
  #                     "APOBEC")]
  
  
  #### . . . Removing the GISTIC peaks Amp_14q11.2 and Del_21p11.2
  clinic10=clinic9[,c("sample",
                      "CNV_chr3.gain","CNV_chr5.gain",   
                      "CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_chr18.gain",  
                      "CNV_Gain_Amp1q","CNV_Amp_17q22","CNV_Amp_2q24.3","CNV_Amp_2q32.1","CNV_Amp_4q13.2","CNV_Amp_6p24.3",  
                      "CNV_Amp_8q24.21","CNV_chr8p.loss","CNV_Del_10p15.3","CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_11q22.1","CNV_Del_12p13.2", 
                      "CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_20p13","CNV_Del_20q13.12","CNV_Del_22q13.32","CNV_Del_2q31.1",  
                      "CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_7p11.2","CNV_Del_7p22.2","CNV_Del_8q24.21","CNV_Del_9p21.3","CNV.SNV_ARID1A",  
                      "CNV.SNV_ARID2","CNV.SNV_ATM","CNV.SNV_BTG1","CNV.SNV_CDKN1B","CNV.SNV_CDKN2C","CNV.SNV_CREBBP","CNV.SNV_CYLD",    
                      "CNV.SNV_DNMT3A","CNV.SNV_EP300","CNV.SNV_TENT5C","CNV.SNV_FUBP1","CNV.SNV_KMT2B","CNV.SNV_KMT2C","CNV.SNV_MAX",     
                      "CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_NFKBIA","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_RASA2","CNV.SNV_RB1",     
                      "CNV.SNV_RPL5","CNV.SNV_SETD2","CNV.SNV_SP140","CNV.SNV_TET2","CNV.SNV_TGDS","CNV.SNV_TP53","CNV.SNV_TRAF2",   
                      "CNV.SNV_TRAF3","CNV.SNV_UBR5","CNV.SNV_ZNF292","SNV_BRAF","SNV_CCND1","SNV_DIS3","SNV_FGFR3",       
                      "SNV_H1.4","SNV_IDH1","SNV_IDH2","SNV_IRF4","SNV_KRAS","SNV_NFKB2","SNV_NRAS",        
                      "SNV_PIK3CA","SNV_PTPN11","SNV_SF3B1","SNV_ABCF1","SNV_ACTG1","SNV_BHLHE41","SNV_DDX3X",       
                      "SNV_DTX1","SNV_DUSP2","SNV_H1.5","SNV_H1.2","SNV_H1.3","SNV_H2BC12","SNV_HLA.C",       
                      "SNV_HUWE1","SNV_INO80E","SNV_IRF1","SNV_KLHL6","SNV_LCE1D","SNV_LTB","SNV_MAN2C1",      
                      "SNV_MRGPRX4","SNV_NDUFC2","SNV_PABPC1","SNV_PANK3","SNV_PRKD2","SNV_RFTN1","SNV_RPRD1B",      
                      "SNV_RPS3A","SNV_SAMHD1","SNV_SLC35G5","SNV_SYT17","SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7",        
                      "SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_BCL7A","SNV_EGR1",        
                      "SNV_MAF","SNV_MAFB","SNV_MAML2","SNV_PIM1","SNV_TCL1A","SNV_XBP1",         
                      "APOBEC")]
  
  write.table(clinic10,
              OUTPUT,
              col.names = T,
              row.names = F,
              quote = F,
              sep = "\t")
  
  # #######
  # #### . . .    - TRANSLOCATION -     . . . ####
  # #######
  # tra=read.delim(paste(INPUT.FILE.FOLDER,"translocations.txt",sep=""),stringsAsFactors = F)
  # rownames(tra)=gsub("_1_BM","",rownames(tra))
  # 
  # setdiff(clinic6$sample,rownames(tra))
  # 
  # # . - Merge CLINIC info with translocation info 
  # clinic7=merge(clinic6,tra,by.x="sample","by.y"="row.names",all.x=T)
  # dim(clinic7)
  # head(clinic7)
  
  # #### . . . - CNVSig - . . . ####
  # cnv.sig=read.delim(CNV.SIG,stringsAsFactors = F)
  # clinic9=merge(clinic9,cnv.sig,by="sample")
}







#### . . . - MatrixGenerator()

#### . . . - Function to create the mutation matrix (Oncogenes, TSGs, Unknown Genes, X genes)
####         matrix = 
dis3_hs=function(MATRIX,CODING.MUT){
  #### . . . - Importing Data - . . . ####
  annot=read.delim(CODING.MUT,
                   stringsAsFactors = F)
  mtx=MATRIX
  DIS3.mut=length(annot[annot$gene=="DIS3","sampleID"])
  if(DIS3.mut>0){
    #### . . . - DIS3 mutations - . . . ####
    ##   . . D479, D488, and R780
    annot_dis<- annot[annot$gene =="DIS3",]
    annot_dis=annot_dis[!(is.na(annot_dis$sampleID)),]
    
    annot_dis$cristal<- 0  
    annot_dis$cristal[grep("D479", annot_dis$aachange)]<- "D479"
    annot_dis$cristal[grep("D488", annot_dis$aachange)]<- "D488"
    annot_dis$cristal[grep("R780", annot_dis$aachange)]<- "R780"
    annot_dis$cristal[annot_dis$impact=="Synonymous"]<-"Synonymous"
    colnames(annot_dis)[1]<-"sample"
    annot_dis2<- annot_dis[,c("sample","cristal", "impact","gene")]
    
    #### . . . - Merging the DIS3 information with the genomic matrix - . . . ####
    int_mtx<- merge(annot_dis2, mtx,by="sample")
    int_mtx$cristal_code<-1
    int_mtx$cristal_code[int_mtx$cristal==0]<- 0
    hs.dis3.samples=unique(int_mtx[int_mtx$cristal_code==1,"sample"])
    
    #### . . . HO & noHS column generation - . . . ####
    mtx$SNV.HS_DIS3=ifelse(mtx$sample %in% hs.dis3.samples,1,
                           ifelse(is.na(mtx$SNV_DIS3),NA,0))
    #table(mtx$SNV_DIS3,mtx$SNV.HS_DIS3)
    
    mtx$SNV.noHS_DIS3=ifelse(mtx$sample %in% hs.dis3.samples,0,mtx$SNV_DIS3)
    #table(mtx$SNV_DIS3,mtx$SNV.noHS_DIS3)
    
    output=mtx[,c("sample","SNV.HS_DIS3","SNV.noHS_DIS3")]
  }else{
    print("No mutations in DIS3 - Creation of SNV.HS_DIS3 and SNV.noHS_DIS3 columns using SNV.DIS3 column")
    output=mtx[,c("sample","SNV_DIS3","SNV_DIS3")]
    colnames(output)=c("sample","SNV.HS_DIS3","SNV.noHS_DIS3")
  }
  return(output)
  
}


#### . . . - MMGenoCla()

#### . . . - Function to create the mutation matrix (Oncogenes, TSGs, Unknown Genes, X genes)
####         matrix = 
MMGenoCla=function(input) {
  #### . . . - Preparation Data - . . . ####
  if(class(input)=="character"){
    all=read.delim(input,
                   stringsAsFactors= F, header=T)
  }else{
    all=input
  }
  
  #### . . . - Column annotation - . . . ####
  cat("Example column annotations:", sep="\n")
  cat("CNV: CNV_chr3.gain, CNV_chr9.gain, CNV_Gain_Amp1q, CNV_Amp_8q24.21, CNV_chr8p.loss, CNV_Del_21p11.2, ...",sep="\n")
  cat("TSG: CNV.SNV_TENT5C, CNV.SNV_RPL5, CNV.SNV_RB1, CNV.SNV_TP53, ...",sep="\n")
  cat("Oncogenes & Other mutations: SNV_NRAS, SNV_BRAF, SNV_KRAS, SNV_KDM6A, ...",sep="\n")
  cat("Other drivers: CNV.Sig, APOBEC, t_CCND1, t_MMSET, t_MAF",sep="\n")
  
  answer=readline(prompt="Do you want to proceed? ")
  if(answer %in% c("Yes","Y","y")){
    
    # if(dis3==F){
    #   
    # }
    
    
    #### . . . - Defining output object - . . . #### 
    output=list()
    
    #### . . . - Change translocation, FAM46C and Chromothripsis names - . . . #### 
    colnames(all)=ifelse(colnames(all)=="t(11;14)","t_CCND1",
                         ifelse(colnames(all) %in% c("t(4;14)","t_NSD2"),"t_MMSET",
                                ifelse(colnames(all) %in% c("t(14;16)","t(14;20)"),"t_MAF",colnames(all))))
    colnames(all)[grep("TENT5C",colnames(all))]=ifelse(colnames(all)[grep("TENT5C",colnames(all))]=="CNV.SNV_TENT5C","CNV.SNV_FAM46C",colnames(all))
    colnames(all)[grep("Chromothripsis",colnames(all))]=ifelse(colnames(all)[grep("Chromothripsis",colnames(all))]=="Chromothripsis","CNV.Sig",colnames(all))
    
    #### . . . - Change translocation values - . . . ####
    translocations=c(grep("t_CCND1",colnames(all)),
                     grep("t_MMSET",colnames(all)),
                     grep("t_MAF",colnames(all)))
    for(i in translocations){
      all[,i]=ifelse(all[,i]==1,2,all[,i])
    }
    
    #### . . . - Only Genomic Features - #### 
    copynumbers=colnames(all)[grep("CNV_",colnames(all))]
    tumorsuppressorgenes=colnames(all)[grep("CNV.SNV_",colnames(all))]
    oncogenesothers=c(colnames(all)[grep("^SNV_",colnames(all))],colnames(all)[grep("^SNV.",colnames(all))])
    all=all[,c("sample",copynumbers,tumorsuppressorgenes,oncogenesothers,"CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF")] # Subset with only genomic features
    rownames(all)=all$sample
    all=all[,-1]
    
    all=as.data.frame(all)
    
    #### . . . - Checking the genomic variables - . . . ####
    # CNV=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_chr18.gain",  
    #       "CNV_Gain_Amp1q","CNV_Amp_14q11.2","CNV_Amp_17q22","CNV_Amp_2q24.3","CNV_Amp_2q32.1","CNV_Amp_4q13.2","CNV_Amp_6p24.3","CNV_Amp_8q24.21","CNV_chr8p.loss",
    #       "CNV_Del_10p15.3","CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_11q22.1","CNV_Del_12p13.2","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_20p13","CNV_Del_20q13.12",
    #       "CNV_Del_21p11.2","CNV_Del_22q13.32","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_7p11.2","CNV_Del_7p22.2","CNV_Del_8q24.21","CNV_Del_9p21.3")
    CNV=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_chr18.gain",  
          "CNV_Gain_Amp1q","CNV_Amp_17q22","CNV_Amp_2q24.3","CNV_Amp_2q32.1","CNV_Amp_4q13.2","CNV_Amp_6p24.3","CNV_Amp_8q24.21","CNV_chr8p.loss",
          "CNV_Del_10p15.3","CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_11q22.1","CNV_Del_12p13.2","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_20p13","CNV_Del_20q13.12",
          "CNV_Del_22q13.32","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_7p11.2","CNV_Del_7p22.2","CNV_Del_8q24.21","CNV_Del_9p21.3")
    TSG=c("CNV.SNV_ARID1A","CNV.SNV_ARID2","CNV.SNV_ATM","CNV.SNV_BTG1","CNV.SNV_CDKN1B","CNV.SNV_CDKN2C","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_EP300", 
          "CNV.SNV_FAM46C","CNV.SNV_FUBP1","CNV.SNV_KMT2B","CNV.SNV_KMT2C","CNV.SNV_MAX","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_NFKBIA","CNV.SNV_POT1","CNV.SNV_PRDM1", 
          "CNV.SNV_RASA2","CNV.SNV_RB1","CNV.SNV_RPL5","CNV.SNV_SETD2","CNV.SNV_SP140","CNV.SNV_TET2","CNV.SNV_TGDS","CNV.SNV_TP53","CNV.SNV_TRAF2","CNV.SNV_TRAF3", 
          "CNV.SNV_UBR5","CNV.SNV_ZNF292")
    SNV=c("SNV_BRAF","SNV_CCND1","SNV_DIS3","SNV_FGFR3","SNV_H1.4","SNV_IDH1","SNV_IDH2","SNV_IRF4","SNV_KRAS","SNV_NFKB2",    
          "SNV_NRAS","SNV_PIK3CA","SNV_PTPN11","SNV_SF3B1","SNV_ABCF1","SNV_ACTG1","SNV_BHLHE41","SNV_DDX3X","SNV_DTX1","SNV_DUSP2",    
          "SNV_H1.5","SNV_H1.2","SNV_H1.3","SNV_H2BC12","SNV_HLA.C","SNV_HUWE1","SNV_INO80E","SNV_IRF1","SNV_KLHL6","SNV_LCE1D",    
          "SNV_LTB","SNV_MAN2C1","SNV_MRGPRX4","SNV_NDUFC2","SNV_PABPC1","SNV_PANK3","SNV_PRKD2","SNV_RFTN1","SNV_RPRD1B","SNV_RPS3A",    
          "SNV_SAMHD1","SNV_SLC35G5","SNV_SYT17","SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7","SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A",    
          "SNV_RPL10","SNV_BCL7A","SNV_EGR1","SNV_MAF","SNV_MAFB","SNV_MAML2","SNV_PIM1","SNV_TCL1A","SNV_XBP1","SNV.HS_DIS3","SNV.noHS_DIS3")
    OTHER=c( "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF")
    
    
    if(length(setdiff(c(CNV,TSG,SNV,OTHER),colnames(all)))!=0){
      cat("Missing genomic features: ",sep="\n")
      print(setdiff(c(CNV,TSG,SNV,OTHER),colnames(all)))
      answer=readline(prompt="Do you want to proceed with the classification? ")
      if(answer %in% c("Yes","Y","y")){
        #### . . . - Remove features - . . . ####
        # . . Based on 1933 samples
        all=all[,!colnames(all) %in% c("SNV_IDH1","SNV_IDH2","SNV_PIK3CA","SNV_ABCF1","SNV_BHLHE41","SNV_DDX3X","SNV_DTX1","SNV_H1.5","SNV_H1.3","SNV_H2BC12","SNV_INO80E",   
                                       "SNV_IRF1","SNV_LCE1D","SNV_MAN2C1","SNV_MRGPRX4","SNV_NDUFC2","SNV_PANK3","SNV_RFTN1","SNV_RPRD1B","SNV_RPS3A","SNV_SLC35G5","SNV_SYT17",    
                                       "SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7","SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_MAF","SNV_MAFB","SNV_MAML2",  
                                       "SNV_PIM1","SNV_TCL1A")]
        
        #### . . . - Creating the matrix with samples with all genomic information (Removing the NAs) - . . . ####
        all=all[complete.cases(all),]
        
        #### . . . - Creation combined features (HRD, gains, deletions, 1q)
        all_gen=all
        all_gen$sample=rownames(all_gen)
        all_gen=all_gen[,c(ncol(all_gen),1:(ncol(all_gen)-1))]
        all_gen$HRD<-0
        all_gen$HRD[all_gen$CNV_chr3.gain+ all_gen$CNV_chr5.gain+all_gen$CNV_chr7.gain+
                      all_gen$CNV_chr9.gain + all_gen$CNV_chr11.gain + all_gen$CNV_chr15.gain +
                      all_gen$CNV_chr19.gain + all_gen$CNV_chr21.gain>=2]<- 1
        
        all_gen$add_gains<-0
        all_gen$add_gains[all_gen$CNV_Amp_14q11.2+ all_gen$CNV_Amp_17q22+all_gen$CNV_Amp_2q24.3+
                            all_gen$CNV_Amp_2q32.1 + all_gen$CNV_Amp_4q13.2 + all_gen$CNV_Amp_6p24.3 +
                            all_gen$CNV_Amp_8q24.21 + all_gen$CNV_chr18.gain>=3]<- 1
        
        all_gen$deletions1=0
        col.del=c(grep("CNV.SNV_",colnames(all_gen)),grep("Del",colnames(all_gen)),grep("loss",colnames(all_gen)))
        minor=min(col.del)
        major=max(col.del)
        new.column=ncol(all_gen)
        for(i in 1:length(all_gen$sample)){
          all_gen[i,new.column]=sum(all_gen[i,c(minor:major)],na.rm = T)
        }
        all_gen$deletions=ifelse(all_gen$deletions1>=5,1,0) # The 10% of the number of the features analyzed (Losses and TSGs). They are 50, so the threshold is 5.
        all_gen=all_gen[,-grep("deletions1",colnames(all_gen))]
        
        #all_gen=as.data.frame.matrix(all_gen)
        all_gen$x1q=all_gen$CNV_Gain_Amp1q
        
        #### . . . - Genomic group creation - . . . ####
        all_gen$clusters<- 0
        all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_KRAS>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
        all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_NRAS>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
        all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_BRAF>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
        all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_FGFR3>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 &  all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
        
        all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==0 & all_gen$add_gains==1 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "B_HRD_gains"
        all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==1 & all_gen$add_gains==0 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "C_HRD_deletions"
        
        all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$add_gains==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"

        #all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==0 & all_gen$add_gains==0 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$CNV.Sig==0 & all_gen$x1q==0]<- "A_HRD_RAS"
        #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$deletions==0 & all_gen$add_gains==1 & all_gen$t_MMSET==0 & all_gen$x1q==0]<- "B_HRD_gains"
        #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$deletions==0 & all_gen$add_gains==0 & all_gen$t_MMSET==0 & all_gen$x1q>0]<- "C_HRD_deletions"
        #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$CNV.Sig==1]<- "C_HRD_deletions"
        #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters %in% c("A_HRD_RAS","B_HRD_gains","0_HRD_Simple") & all_gen$x1q>=1 ]<-"C_HRD_deletions"
        #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters %in% c("A_HRD_RAS","B_HRD_gains","0_HRD_Simple") & all_gen$CNV.SNV_TP53>1 ]<-"C_HRD_deletions"
        
        
        
        
        #table(all_gen$clusters, all_gen$HRD)
        all_gen$clusters[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC>1]<-"I_MAF_APOBEC"
        all_gen$clusters[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC<2]<-"I_MAF_APOBEC"
        
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$HRD==1 & all_gen$deletions==0]<-"E_CCND1_HRD"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_KRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_BRAF>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_NRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_FGFR3>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$deletions>0 ]<-"G_CCND1_complex"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters==0 & all_gen$deletions==0]<-"FF_CCND1_RAS"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$x1q>=1 ]<-"G_CCND1_complex"
        all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$CNV.SNV_TP53>1 ]<-"G_CCND1_complex"
        
        #table(all_gen$clusters)
        
        all_gen$clusters[all_gen$t_MMSET>0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"N_NSD2_1q_del13q"
        all_gen$clusters[all_gen$t_MMSET>0 & all_gen$x1q==0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"P_NSD2_del13q"
        all_gen$clusters[all_gen$t_MMSET==0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0) & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"Q_x1q_del13q"
        
        all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_MAF>0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
        all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_MMSET==0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
        all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_CCND1==0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
        
        all_gen$clusters[all_gen$APOBEC>1 & all_gen$clusters == "Q_x1q_del13q"]<-"L_MAF_APOBEC"
        all_gen$clusters[all_gen$APOBEC>1 & 
                           (all_gen$clusters == "A_HRD_RAS" | all_gen$clusters %in% c("B_HRD_gains","BB_HRD_gains") | 
                              all_gen$clusters %in% c("C_HRD_deletions","CC_HRD_deletions") | all_gen$clusters %in% c("D_HRD_other","DD_HRD_other") | 
                              all_gen$clusters %in% c("0_HRD_Simple"))]<-"M_MAF_APOBEC"
        
        all_gen$clusters[all_gen$t_MMSET>0 & all_gen$clusters == 0]<-"O_NSD2"
        
        #table(all_gen[all_gen$clusters==0,"deletions"])
        all_gen$clusters[all_gen$clusters==0 & all_gen$deletions==1]<-"R_Multiple_losses"
        #all_gen[all_gen$clusters==0,]
        all_gen$clusters[all_gen$clusters==0]<-"RR_Simple"
        #table(all_gen$clusters)
        
        all_gen$clusters[all_gen$clusters==0 & all_gen$t_CCND1>0]<-"HH_CCND1"
        
        all_gen$clusters[all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") &
                           (all_gen$APOBEC>1 | all_gen$CNV.Sig==1)]<-"G_CCND1_complex"
        
        all_gen$clusters[all_gen$clusters %in% c("M_MAF_APOBEC","I_MAF_APOBEC","L_MAF_APOBEC") &
                           (all_gen$t_MMSET>0 & all_gen$t_MAF>0)]<-"N_NSD2_1q_del13q"
        
        all_gen$clusters[all_gen$clusters %in% c("RR_Simple") &
                           (all_gen$CNV.Sig==1 | all_gen$x1q==2)]<-"R_Multiple_losses"
        
        # all_gen$clusters[all_gen$clusters %in% c("A_HRD_RAS","B_HRD_gains") &
        #                    (all_gen$Gain_Amp1q==2)]<-"C_HRD_deletions"
        
        #### . . . - Changing clomun names for the heatmap - . . . ####
        all_gen<- all_gen[order(all_gen$clusters),]
        
        colnames(all_gen)=gsub("CNV_","",colnames(all_gen))
        colnames(all_gen)=gsub("CNV.SNV_","",colnames(all_gen))
        colnames(all_gen)=gsub("SNV_","",colnames(all_gen))
        colnames(all_gen)=gsub("SNV.","",colnames(all_gen))
        
        colnames(all_gen)[grep("t_CCND1",colnames(all_gen))]="CCND1-IGH"
        colnames(all_gen)[grep("t_MMSET",colnames(all_gen))]="NSD2-IGH"
        colnames(all_gen)[grep("t_MAF",colnames(all_gen))]="MAF-IGH"
        
        #### . . . - Sorting the matrix for the heatmap - . . . ####
        all_gen1=all_gen[,c("sample", #samples
                            "chr3.gain","chr5.gain","chr7.gain","chr9.gain","chr11.gain","chr15.gain","chr19.gain","chr21.gain", #HRD
                            "Gain_Amp1q", #1q
                            "chr18.gain",#chr18.gain +
                            #"Amp_14q11.2",
                            "Amp_17q22","Amp_2q24.3","Amp_2q32.1","Amp_4q13.2","Amp_6p24.3","Amp_8q24.21",  # other gains
                            "ARID1A","CDKN2C","FAM46C","FUBP1","RPL5", #chr1
                            "Del_2q31.1","Del_2q37.3","DNMT3A","SP140", #chr2
                            "RASA2","SETD2", #chr3
                            "TET2", #chr4
                            "Del_6q26","PRDM1","ZNF292", #chr6
                            "Del_7p11.2","Del_7p22.2","KMT2C","POT1", #chr7
                            "chr8p.loss","Del_8q24.21","UBR5", # chr8
                            "Del_9p21.3","TRAF2", #chr9
                            "Del_10p15.3","Del_10q24.32","Del_10q26.3", #chr10
                            "Del_11q22.1","ATM", #chr11
                            "Del_12p13.2","Del_12q24.31","ARID2","BTG1","CDKN1B", #chr12
                            "RB1","TGDS", #chr13
                            "Del_14q24.3","MAX","NFKBIA","TRAF3", #chr14
                            "CREBBP","CYLD", #chr16
                            "NCOR1","NF1","TP53", #chr17
                            "KMT2B", #chr19
                            "Del_20p13","Del_20q13.12", #chr20
                            #"Del_21p11.2", #chr21
                            "Del_22q13.32","EP300", #chr22
                            "BRAF","KRAS","NRAS", #BRAF,KRAS,NRAS # (MAPK pathway)
                            "CCND1","FGFR3","H1.4","IRF4","NFKB2","PTPN11","SF3B1","ACTG1","DUSP2","H1.2","HLA.C",
                            "HUWE1","KLHL6","LTB","PABPC1","PRKD2","SAMHD1","BCL7A","EGR1","XBP1",#other mutations  # N.B.: Remove DIS3
                            "HS_DIS3","noHS_DIS3", # DIS3 HS,DIS3 noHS
                            "CCND1-IGH","NSD2-IGH","MAF-IGH", #translocations
                            "APOBEC", #APOBEC
                            "CNV.Sig", #CNV.Sig
                            "HRD","clusters"#"TC6_calls","RNA_high_risk","UAMS.call" [These classification are not present in this code]
        )]
        
        annotation_col_genomic<- as.data.frame((all_gen1[,c("clusters")]))
        rownames(annotation_col_genomic)<-all_gen1$sample
        colnames(annotation_col_genomic)<-c("clusters")
        
        #### . . . - Defining new cluster names -  . . . ####
        annotation_col_genomic$clusters.new=ifelse(annotation_col_genomic$clusters %in% c("A_HRD_RAS"),"HRD_RAS",
                                                   ifelse(annotation_col_genomic$clusters %in% c("B_HRD_gains"),"HRD_Gains",
                                                          ifelse(annotation_col_genomic$clusters %in% c("C_HRD_deletions"),"HRD_Complex_Cytogenetic",
                                                                 ifelse(annotation_col_genomic$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS"),"CCND1_Simple",
                                                                        ifelse(annotation_col_genomic$clusters %in% c("G_CCND1_complex"),"CCND1_Complex_Cytogenetic",
                                                                               ifelse(annotation_col_genomic$clusters %in% c("I_MAF_APOBEC","L_MAF_APOBEC","M_MAF_APOBEC"),"MAF_and_or_HyperAPOBEC",
                                                                                      ifelse(annotation_col_genomic$clusters %in% c("N_NSD2_1q_del13q"),"NSD2_GainAmp1q_Del13q",
                                                                                             ifelse(annotation_col_genomic$clusters %in% c("O_NSD2"),"NSD2_HRD",
                                                                                                    ifelse(annotation_col_genomic$clusters %in% c("P_NSD2_del13q"),"NSD2_Del13q",
                                                                                                           ifelse(annotation_col_genomic$clusters %in% c("Q_x1q_del13q"),"GainAmp1q_Del13q",
                                                                                                                  ifelse(annotation_col_genomic$clusters %in% c("R_Multiple_losses"),"Multiple_Losses",
                                                                                                                         ifelse(annotation_col_genomic$clusters %in% c("RR_Simple"),"Simple",NA))))))))))))
        
        
        mycol_plus<- c(brewer.pal(8,"Dark2"),brewer.pal(8,"Paired"), brewer.pal(8,"Set2"))
        ann_colors = list(
          clusters.new =c("HRD_RAS"="lightpink3",
                          "HRD_Gains"="mediumpurple2",
                          "HRD_Complex_Cytogenetic"="darkorchid4",
                          "CCND1_Simple"="skyblue3",
                          "CCND1_Complex_Cytogenetic"="navy",
                          "MAF_and_or_HyperAPOBEC"="darkgreen",
                          "NSD2_GainAmp1q_Del13q"="lightgoldenrod",
                          "NSD2_HRD"="gold",
                          "NSD2_Del13q"="goldenrod4",
                          "GainAmp1q_Del13q"="sienna3",
                          "Multiple_Losses"="firebrick",
                          "Simple"="slategray"
          )
        )
        
        
        # space_heat<- as.numeric(table(all_gen1$clusters))
        # 
        # space_heat2<- c(space_heat[1],
        #                 space_heat[1] + space_heat[2],
        #                 space_heat[1] + space_heat[2] + space_heat[3],
        #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4],
        #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5],
        #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]),
        #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8],
        #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8] + space_heat[9]
        # )
        
        # space.row=c(8,
        #             8 + 1,
        #             8 + 1 + 8,
        #             8 + 1 + 8 + 51,
        #             8 + 1 + 8 + 51 + 3,
        #             8 + 1 + 8 + 51 + 3 + 22,
        #             8 + 1 + 8 + 51 + 3 + 22 + 3,
        #             8 + 1 + 8 + 51 + 3 + 22 + 3 + 1,
        #             8 + 1 + 8 + 51 + 3 + 22 + 3 + 1 + 1)
        
        annotation_col_genomic1=as.data.frame(annotation_col_genomic[,c("clusters.new")])
        colnames(annotation_col_genomic1)="clusters.new"
        rownames(annotation_col_genomic1)=rownames(annotation_col_genomic)
        rownames(all_gen1)=all_gen1$sample
        all_gen1=all_gen1[,-grep("sample",colnames(all_gen1))]
        mtx=all_gen1[,-c(grep("HRD",colnames(all_gen1)),grep("clusters",colnames(all_gen1)))]
        
        output$genomic=mtx
        
        annotation_col_genomic=annotation_col_genomic[order(annotation_col_genomic$clusters),]
        output$annotation=annotation_col_genomic1
        
        mtx=mtx[match(rownames(annotation_col_genomic),rownames(mtx)),]
        #identical(rownames(mtx),rownames(annotation_col_genomic))
        
        output$colors=ann_colors
        return(output)
        
      }else{
        cat("STOP. Check the genomic features before to run the genomic classification",sep="\n")
      }
    }else{
      answer=readline(prompt="The matrix is ready for the classification, do you want to proceed? ")
    }
    if(answer %in% c("Yes","Y","y")){
      
      #### . . . - Remove features - . . . ####
      # . . Based on 1933 samples
      high.missing.features=all[,colnames(all) %in% c("sample","SNV_IDH1","SNV_IDH2","SNV_PIK3CA","SNV_ABCF1","SNV_BHLHE41","SNV_DDX3X","SNV_DTX1","SNV_H1.5","SNV_H1.3","SNV_H2BC12","SNV_INO80E",   
                                                     "SNV_IRF1","SNV_LCE1D","SNV_MAN2C1","SNV_MRGPRX4","SNV_NDUFC2","SNV_PANK3","SNV_RFTN1","SNV_RPRD1B","SNV_RPS3A","SNV_SLC35G5","SNV_SYT17",    
                                                     "SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7","SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_MAF","SNV_MAFB","SNV_MAML2",  
                                                     "SNV_PIM1","SNV_TCL1A")]
      
      all=all[,!colnames(all) %in% c("SNV_IDH1","SNV_IDH2","SNV_PIK3CA","SNV_ABCF1","SNV_BHLHE41","SNV_DDX3X","SNV_DTX1","SNV_H1.5","SNV_H1.3","SNV_H2BC12","SNV_INO80E",   
                                     "SNV_IRF1","SNV_LCE1D","SNV_MAN2C1","SNV_MRGPRX4","SNV_NDUFC2","SNV_PANK3","SNV_RFTN1","SNV_RPRD1B","SNV_RPS3A","SNV_SLC35G5","SNV_SYT17",    
                                     "SNV_TBC1D29P","SNV_TNFRSF11B","SNV_USP7","SNV_ZFP36L1","SNV_ATRX","SNV_KDM5C","SNV_KDM6A","SNV_RPL10","SNV_MAF","SNV_MAFB","SNV_MAML2",  
                                     "SNV_PIM1","SNV_TCL1A")]
      
      #### . . . - Creating the matrix with samples with all genomic information (Removing the NAs) - . . . ####
      all=all[complete.cases(all),]
      
      #### . . . - Creation combined features (HRD, gains, deletions, 1q)
      all_gen=all
      all_gen$sample=rownames(all_gen)
      all_gen=all_gen[,c(ncol(all_gen),1:(ncol(all_gen)-1))]
      all_gen$HRD<-0
      all_gen$HRD[all_gen$CNV_chr3.gain+ all_gen$CNV_chr5.gain+all_gen$CNV_chr7.gain+
                    all_gen$CNV_chr9.gain + all_gen$CNV_chr11.gain + all_gen$CNV_chr15.gain +
                    all_gen$CNV_chr19.gain + all_gen$CNV_chr21.gain>=2]<- 1
      
      all_gen$add_gains<-0
      all_gen$add_gains[all_gen$CNV_Amp_14q11.2+ all_gen$CNV_Amp_17q22+all_gen$CNV_Amp_2q24.3+
                          all_gen$CNV_Amp_2q32.1 + all_gen$CNV_Amp_4q13.2 + all_gen$CNV_Amp_6p24.3 +
                          all_gen$CNV_Amp_8q24.21 + all_gen$CNV_chr18.gain>=3]<- 1
      
      all_gen$deletions1=0
      col.del=c(grep("CNV.SNV_",colnames(all_gen)),grep("Del",colnames(all_gen)),grep("loss",colnames(all_gen)))
      minor=min(col.del)
      major=max(col.del)
      new.column=ncol(all_gen)
      for(i in 1:length(all_gen$sample)){
        all_gen[i,new.column]=sum(all_gen[i,c(minor:major)],na.rm = T)
      }
      all_gen$deletions=ifelse(all_gen$deletions1>=5,1,0) # The 10% of the number of the features analyzed (Losses and TSGs). They are 50, so the threshold is 5.
      all_gen=all_gen[,-grep("deletions1",colnames(all_gen))]
      
      #all_gen=as.data.frame.matrix(all_gen)
      all_gen$x1q=all_gen$CNV_Gain_Amp1q
      
      #### . . . - Genomic group creation - . . . ####
      all_gen$clusters<- 0
      all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_KRAS>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
      all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_NRAS>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
      all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_BRAF>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
      all_gen$clusters[all_gen$HRD==1 & all_gen$SNV_FGFR3>0 & all_gen$add_gains==0 & all_gen$deletions==0 & all_gen$t_MMSET==0 &  all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"
      
      all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==0 & all_gen$add_gains==1 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "B_HRD_gains"
      all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==1 & all_gen$add_gains==0 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "C_HRD_deletions"
      
      all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$add_gains==0 & all_gen$t_MMSET==0 & all_gen$t_CCND1==0 & all_gen$t_MAF==0]<- "A_HRD_RAS"

      #all_gen$clusters[all_gen$HRD==1 & all_gen$deletions==0 & all_gen$add_gains==0 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$CNV.Sig==0 & all_gen$x1q==0]<- "A_HRD_RAS"
      #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$deletions==0 & all_gen$add_gains==1 & all_gen$t_MMSET==0 & all_gen$x1q==0]<- "B_HRD_gains"
      #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$deletions==0 & all_gen$add_gains==0 & all_gen$t_MMSET==0 & all_gen$x1q>0]<- "C_HRD_deletions"
      #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters==0 & all_gen$t_MMSET==0 & all_gen$CNV.Sig==1]<- "C_HRD_deletions"
      #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters %in% c("A_HRD_RAS","B_HRD_gains","0_HRD_Simple") & all_gen$x1q>=1 ]<-"C_HRD_deletions"
      #all_gen$clusters[all_gen$HRD==1 & all_gen$clusters %in% c("A_HRD_RAS","B_HRD_gains","0_HRD_Simple") & all_gen$CNV.SNV_TP53>1 ]<-"C_HRD_deletions"
      
      
      
      
      #table(all_gen$clusters, all_gen$HRD)
      all_gen$clusters[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC>1]<-"I_MAF_APOBEC"
      all_gen$clusters[all_gen$t_MAF>0 & all_gen$SNV.HS_DIS3==0 & all_gen$APOBEC<2]<-"I_MAF_APOBEC"
      
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$HRD==1 & all_gen$deletions==0]<-"E_CCND1_HRD"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_KRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_BRAF>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_NRAS>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$SNV_FGFR3>0 & all_gen$deletions==0]<-"F_CCND1_RAS"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$deletions>0 ]<-"G_CCND1_complex"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters==0 & all_gen$deletions==0]<-"FF_CCND1_RAS"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$x1q>=1 ]<-"G_CCND1_complex"
      all_gen$clusters[all_gen$t_CCND1>0 & all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") & all_gen$CNV.SNV_TP53>1 ]<-"G_CCND1_complex"
      
      #table(all_gen$clusters)
      
      all_gen$clusters[all_gen$t_MMSET>0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"N_NSD2_1q_del13q"
      all_gen$clusters[all_gen$t_MMSET>0 & all_gen$x1q==0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0)  & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"P_NSD2_del13q"
      all_gen$clusters[all_gen$t_MMSET==0 & all_gen$x1q>0 & (all_gen$CNV.SNV_RB1>0 | all_gen$CNV.SNV_TGDS>0) & all_gen$SNV.HS_DIS3==0 & all_gen$clusters==0]<-"Q_x1q_del13q"
      
      all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_MAF>0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
      all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_MMSET==0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
      all_gen$clusters[all_gen$APOBEC>1 & all_gen$t_CCND1==0 & all_gen$clusters == 0]<-"L_MAF_APOBEC"
      
      all_gen$clusters[all_gen$APOBEC>1 & all_gen$clusters == "Q_x1q_del13q"]<-"L_MAF_APOBEC"
      all_gen$clusters[all_gen$APOBEC>1 & 
                         (all_gen$clusters == "A_HRD_RAS" | all_gen$clusters %in% c("B_HRD_gains","BB_HRD_gains") | all_gen$clusters %in% c("C_HRD_deletions","CC_HRD_deletions") | all_gen$clusters %in% c("D_HRD_other","DD_HRD_other") | all_gen$clusters %in% c("0_HRD_Simple"))]<-"M_MAF_APOBEC"
      
      all_gen$clusters[all_gen$t_MMSET>0 & all_gen$clusters == 0]<-"O_NSD2"
      
      #table(all_gen[all_gen$clusters==0,"deletions"])
      all_gen$clusters[all_gen$clusters==0 & all_gen$deletions==1]<-"R_Multiple_losses"
      #all_gen[all_gen$clusters==0,]
      all_gen$clusters[all_gen$clusters==0]<-"RR_Simple"
      #table(all_gen$clusters)
      
      all_gen$clusters[all_gen$clusters==0 & all_gen$t_CCND1>0]<-"HH_CCND1"
      
      all_gen$clusters[all_gen$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS") &
                         (all_gen$APOBEC>1 | all_gen$CNV.Sig==1)]<-"G_CCND1_complex"
      
      all_gen$clusters[all_gen$clusters %in% c("M_MAF_APOBEC","I_MAF_APOBEC","L_MAF_APOBEC") &
                         (all_gen$t_MMSET>0 & all_gen$t_MAF>0)]<-"N_NSD2_1q_del13q"
      
      all_gen$clusters[all_gen$clusters %in% c("RR_Simple") &
                         (all_gen$CNV.Sig==1 | all_gen$x1q==2)]<-"R_Multiple_losses"
      
      
      #### . . . - Changing clomun names for the heatmap - . . . ####
      all_gen<- all_gen[order(all_gen$clusters),]
      
      colnames(all_gen)=gsub("CNV_","",colnames(all_gen))
      colnames(all_gen)=gsub("CNV.SNV_","",colnames(all_gen))
      colnames(all_gen)=gsub("SNV_","",colnames(all_gen))
      colnames(all_gen)=gsub("SNV.","",colnames(all_gen))
      
      colnames(all_gen)[grep("t_CCND1",colnames(all_gen))]="CCND1-IGH"
      colnames(all_gen)[grep("t_MMSET",colnames(all_gen))]="NSD2-IGH"
      colnames(all_gen)[grep("t_MAF",colnames(all_gen))]="MAF-IGH"
      
      #### . . . - Sorting the matrix for the heatmap - . . . ####
      all_gen1=all_gen[,c("sample", #samples
                          "chr3.gain","chr5.gain","chr7.gain","chr9.gain","chr11.gain","chr15.gain","chr19.gain","chr21.gain", #HRD
                          "Gain_Amp1q", #1q
                          "chr18.gain", #chr18.gain + 
                          #"Amp_14q11.2",
                          "Amp_17q22","Amp_2q24.3","Amp_2q32.1","Amp_4q13.2","Amp_6p24.3","Amp_8q24.21", #other gains
                          "ARID1A","CDKN2C","FAM46C","FUBP1","RPL5", #chr1
                          "Del_2q31.1","Del_2q37.3","DNMT3A","SP140", #chr2
                          "RASA2","SETD2", #chr3
                          "TET2", #chr4
                          "Del_6q26","PRDM1","ZNF292", #chr6
                          "Del_7p11.2","Del_7p22.2","KMT2C","POT1", #chr7
                          "chr8p.loss","Del_8q24.21","UBR5", # chr8
                          "Del_9p21.3","TRAF2", #chr9
                          "Del_10p15.3","Del_10q24.32","Del_10q26.3", #chr10
                          "Del_11q22.1","ATM", #chr11
                          "Del_12p13.2","Del_12q24.31","ARID2","BTG1","CDKN1B", #chr12
                          "RB1","TGDS", #chr13
                          "Del_14q24.3","MAX","NFKBIA","TRAF3", #chr14
                          "CREBBP","CYLD", #chr16
                          "NCOR1","NF1","TP53", #chr17
                          "KMT2B", #chr19
                          "Del_20p13","Del_20q13.12", #chr20
                          #"Del_21p11.2", #chr21
                          "Del_22q13.32","EP300", #chr22
                          "BRAF","KRAS","NRAS", #BRAF,KRAS,NRAS # (MAPK pathway)
                          "CCND1","FGFR3","H1.4","IRF4","NFKB2","PTPN11","SF3B1","ACTG1","DUSP2","H1.2","HLA.C",
                          "HUWE1","KLHL6","LTB","PABPC1","PRKD2","SAMHD1","BCL7A","EGR1","XBP1",#other mutations  # N.B.: Remove DIS3
                          "HS_DIS3","noHS_DIS3", # DIS3 HS,DIS3 noHS
                          "CCND1-IGH","NSD2-IGH","MAF-IGH", #translocations
                          "APOBEC", #APOBEC
                          "CNV.Sig", #CNV.Sig
                          "HRD","clusters"#"TC6_calls","RNA_high_risk","UAMS.call" [These classification are not present in this code]
      )]
      
      annotation_col_genomic<- as.data.frame((all_gen1[,c("clusters")]))
      rownames(annotation_col_genomic)<-all_gen1$sample
      colnames(annotation_col_genomic)<-c("clusters")
      
      #### . . . - Defining new cluster names -  . . . ####
      annotation_col_genomic$clusters.new=ifelse(annotation_col_genomic$clusters %in% c("A_HRD_RAS"),"HRD_RAS",
                                                 ifelse(annotation_col_genomic$clusters %in% c("B_HRD_gains"),"HRD_Gains",
                                                        ifelse(annotation_col_genomic$clusters %in% c("C_HRD_deletions"),"HRD_Complex_Cytogenetic",
                                                               ifelse(annotation_col_genomic$clusters %in% c("E_CCND1_HRD","F_CCND1_RAS","FF_CCND1_RAS"),"CCND1_Simple",
                                                                      ifelse(annotation_col_genomic$clusters %in% c("G_CCND1_complex"),"CCND1_Complex_Cytogenetic",
                                                                             ifelse(annotation_col_genomic$clusters %in% c("I_MAF_APOBEC","L_MAF_APOBEC","M_MAF_APOBEC"),"MAF_and_or_HyperAPOBEC",
                                                                                    ifelse(annotation_col_genomic$clusters %in% c("N_NSD2_1q_del13q"),"NSD2_GainAmp1q_Del13q",
                                                                                           ifelse(annotation_col_genomic$clusters %in% c("O_NSD2"),"NSD2_HRD",
                                                                                                  ifelse(annotation_col_genomic$clusters %in% c("P_NSD2_del13q"),"NSD2_Del13q",
                                                                                                         ifelse(annotation_col_genomic$clusters %in% c("Q_x1q_del13q"),"GainAmp1q_Del13q",
                                                                                                                ifelse(annotation_col_genomic$clusters %in% c("R_Multiple_losses"),"Multiple_Losses",
                                                                                                                       ifelse(annotation_col_genomic$clusters %in% c("RR_Simple"),"Simple",NA))))))))))))
      
      
      mycol_plus<- c(brewer.pal(8,"Dark2"),brewer.pal(8,"Paired"), brewer.pal(8,"Set2"))
      ann_colors = list(
        clusters.new =c("HRD_RAS"="lightpink3",
                        "HRD_Gains"="mediumpurple2",
                        "HRD_Complex_Cytogenetic"="darkorchid4",
                        "CCND1_Simple"="skyblue3",
                        "CCND1_Complex_Cytogenetic"="navy",
                        "MAF_and_or_HyperAPOBEC"="darkgreen",
                        "NSD2_GainAmp1q_Del13q"="lightgoldenrod",
                        "NSD2_HRD"="gold",
                        "NSD2_Del13q"="goldenrod4",
                        "GainAmp1q_Del13q"="sienna3",
                        "Multiple_Losses"="firebrick",
                        "Simple"="slategray"
        )
      )
      
      
      # space_heat<- as.numeric(table(all_gen1$clusters))
      # 
      # space_heat2<- c(space_heat[1],
      #                 space_heat[1] + space_heat[2],
      #                 space_heat[1] + space_heat[2] + space_heat[3],
      #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4],
      #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5],
      #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]),
      #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8],
      #                 space_heat[1] + space_heat[2] + space_heat[3] + space_heat[4] + space_heat[5] + c(space_heat[6] + space_heat[7]) + space_heat[8] + space_heat[9]
      # )
      
      # space.row=c(8,
      #             8 + 1,
      #             8 + 1 + 8,
      #             8 + 1 + 8 + 51,
      #             8 + 1 + 8 + 51 + 3,
      #             8 + 1 + 8 + 51 + 3 + 22,
      #             8 + 1 + 8 + 51 + 3 + 22 + 3,
      #             8 + 1 + 8 + 51 + 3 + 22 + 3 + 1,
      #             8 + 1 + 8 + 51 + 3 + 22 + 3 + 1 + 1)
      
      annotation_col_genomic1=as.data.frame(annotation_col_genomic[,c("clusters.new")])
      colnames(annotation_col_genomic1)="clusters.new"
      rownames(annotation_col_genomic1)=rownames(annotation_col_genomic)
      rownames(all_gen1)=all_gen1$sample
      all_gen1=all_gen1[,-grep("sample",colnames(all_gen1))]
      mtx=all_gen1[,-c(grep("HRD",colnames(all_gen1)),grep("clusters",colnames(all_gen1)))]
      
      output$genomic=mtx
      
      annotation_col_genomic=annotation_col_genomic[order(annotation_col_genomic$clusters),]
      output$annotation=annotation_col_genomic1
      
      mtx=mtx[match(rownames(annotation_col_genomic),rownames(mtx)),]
      #identical(rownames(mtx),rownames(annotation_col_genomic))
      
      output$colors=ann_colors
      
      output$high.missing.features=high.missing.features
      return(output)
      
    }else{
      cat("Are you sure? The matrix is ready for the genomic classification!",sep="\n")
    } 
  }else{
    cat("STOP",sep="\n")
  }
}

















































