ann_exp_cpm1<-read.csv("ann_exp_cpm1.csv")
ann_met1<-read.csv("ann_met1.csv")
snp_pos<-read.csv("snp_pos.csv")
met_pos<-read.csv("met_pos.csv")
p_d_fin<-read.csv("p_d_fin.csv")


library(Biobase)
library(tidyverse)
library(MatrixEQTL)

ann_exp_cpm1[ann_exp_cpm1 == "AA"]<-0
ann_exp_cpm1[ann_exp_cpm1 == "AB"]<-1
ann_exp_cpm1[ann_exp_cpm1 == "BB"]<-2
ann_exp_cpm1[ann_exp_cpm1 == "NC"]<-NA

head(ann_exp_cpm1)
head(ann_exp_cpm1)
ins_c<-c("rs10176669","rs4438452","rs12446308","rs13041757","rs10517215","rs10895256",
         "rs1820453","rs716274","rs1656402","rs1209950","rs9981861","rs1878022",
         "rs10937823","rs7629386","rs969088","rs41997","rs12000445","rs3850370","rs1571228",
         "rs2371030","rs11098063","rs11568927","rs6901416","rs10766739","rs10023113","rs2107561",
         "rs6882451","rs1826692","rs6595026")
in_s_pos<-snp_pos[snp_pos$SNP_ID %in% ins_c,]

fil1<- ann_exp_cpm1[ann_exp_cpm1$SNP_ID %in% ins_c, ]
head(in_s_pos)
in_s_pos<-in_s_pos
fil1<-fil1%>% select(colnames(fil1)[2:156])

fil1[fil1 == "AA"]<-0
fil1[fil1 == "AB"]<-1
fil1[fil1 == "BB"]<-2
fil1[fil1 == "NC"]<-NA


p_d_fin[p_d_fin == "Male"]<-0.0
p_d_fin[p_d_fin == "Female"]<-1.0
p_d_fin[p_d_fin == "Caucasian"]<-0.0
p_d_fin[p_d_fin == "Asian"]<-1.0
p_d_fin[p_d_fin == "	
American Indian or Alaska Native"]<-2.0

write.table(fil1,"fil11.csv",sep = "\t")
fil1<-read.csv("fil1.csv",header = TRUE)
rownames(fil1)<-fil1$SNP_ID
fil1<-fil1%>% select(-c("SNP_ID"))
rownames(p_d_fin)<-p_d_fin$ID
p_d_fin<-p_d_fin%>% select(-c("ID"))
p_d_fin<-read.table("p_d_fin.csv")
write.table(p_d_fin,"p_d_fin_alt.csv",sep = "\t")
colnames(p_d_fin)[1]<-"ID"

snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
snps$LoadFile("fil11.csv");

snps
cvrt <- SlicedData$new()  
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile("p_d_fin_alt.csv")



cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length("p_d_fin_alt.csv")>0) {
  cvrt$LoadFile("p_d_fin_alt.csv");
}

ann_met1<-ann_met1%>% select(-c("ID"))
rownames(ann_met1)<-ann_met1$ID
write.table(ann_met2,"ann_met2.csv",sep = "\t")


meth <- SlicedData$new()
meth$fileDelimiter <- "\t"
meth$fileOmitCharacters <- "NA"
meth$fileSkipRows <- 1
meth$fileSkipColumns <- 1
meth$fileSliceSize <- 20000
meth$LoadFile("ann_met2.csv")

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

cisDist = 1e6;
useModel = modelLINEAR


me = Matrix_eQTL_main(
  snps = snps,
  gene = meth,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = in_s_pos,
  genepos = met_pos_nar,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
is.numeric(in_s_pos[[3]])
in_s_pos[[3]]
in_s_pos$Physical.Position<-as.numeric(in_s_pos$Physical.Position)
met_pos$RANGE_START<-as.numeric(met_pos$RANGE_START)
met_pos$RANGE_END<-as.numeric(met_pos$RANGE_END)
h<-which(is.na(met_pos$RANGE_START))
length(h)
id_na<-met_pos$ID[which(is.na(met_pos$RANGE_START))]
id_na[1:5]
met_pos_nar<-met_pos[-h,]
ann_met2<-ann_met1[-h,]
  
unlink(output_file_name_tra);
unlink(output_file_name_cis);
out_cis<-as.data.frame(show(me$cis$eqtls))
show(me$trans$eqtls)  
plot(me)  
library(dplyr)
which(ann_exp_cpm1$SNP_ID=="rs13041757")
trans1<-ann_exp_cpm1[283584,]
trans1_cpg<-ann_met1[rownames(ann_met1)=="cg12128876",]

library(tidyverse)
##############################
# wide to long format for plotting the methylation values 
trans1_cpg<-trans1_cpg %>%
  pivot_longer(cols = colnames(trans1_cpg),names_to = "probe",values_to = "met" )
ggplot(a_g,aes(probe,met))+geom_boxplot()
##############################

trans1<-trans1%>%select(-c("X"))
trans1<-trans1 %>%
  pivot_longer(cols = colnames(trans1),names_to = "sample",values_to = "genotype" )

colnames(trans1_cpg)[1]<-"sample"

plot_f<-merge.data.frame(trans1,trans1_cpg,by = "sample")
g1<-ggplot(plot_f,aes(genotype,met),)+geom_boxplot( colour = "black",lwd = 1, flatten = 4,fill = "#c0c0c0")+
  
  geom_jitter(shape=21,aes(fill = genotype),show.legend = T,size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text", x = 3, y =0.1 , label = paste(	
   "rs13041757"  ,"cg12128876",sep = "~"))+xlab("Genotype")+
  ylab("Beta value")
g1

ggplot(plot_f,aes(genotype,met))+geom_boxplot( colour = "black",lwd = 1, flatten = 4,fill = "#e9ffdb")+
  
  geom_jitter(shape=21,aes(fill = genotype),show.legend = F,size = 4)
  



x= "rs1209950"
y ="cg11890956"
ann_exp_cpm1[ann_exp_cpm1 == "NC"]<-NA
ann_exp_cpm1<-ann_exp_cpm1%>%select(-c("X"))
rownames(ann_exp_cpm1)<-ann_exp_cpm1$SNP_ID
?stat_boxplot
?geom_boxplot





########################
plot<-function(a,b){
  v1<-which(ann_exp_cpm1$SNP_ID == a)
  v2<-which(rownames(ann_met1) == b)
  trans_a<-ann_exp_cpm1[v1,]
  trans_b<-ann_met1[v2,]
  trans_a<-trans_a %>%
      pivot_longer(cols = colnames(trans_a),names_to = "sample",values_to = "genotype" )
  trans_a<-trans_a[-c(which(is.na(trans_a$genotype))),]
  trans_b<-trans_b %>%
    pivot_longer(cols = colnames(trans_b),names_to = "sample",values_to = "met" )
  plot_k<-merge.data.frame(trans_a,trans_b,by = "sample")
  g1<-ggplot(plot_k,aes(genotype,met),)+geom_boxplot( colour = "black",lwd = 1, flatten = 4,fill = "#c0c0c0")+
    
    geom_jitter(shape=21,aes(fill = genotype),show.legend = T,size = 4)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    annotate("text", x = 3, y =0.3 , label = paste(a,b,sep = "~"))+xlab("Genotype")+
    ylab("Beta value")

return(g1)
}
################################


plot(	
  "rs17181550",
  "cg03184439")

ii<-plot("rs1878022",
"cg01137198")
itt1<-function(k){
a<- plot(as.character(out_trans_1$snps[k]),
     as.character(out_trans_1$gene[k]))
show(a)

}
itt1(27)
view(ann_exp_cpm1[283584,])
plot("rs17391694",
     "cg24762437")
out_trn<-out_trans_1[ -which(out_trans_1$snps=="rs12000445"),]
out_trn<-out_trn[ -which(out_trn$snps=="rs11568927"),]




















