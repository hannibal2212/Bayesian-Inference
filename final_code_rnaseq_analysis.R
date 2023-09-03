library(dplyr)
library(readr)
library(HTSFilter)
library(biomaRt)
library(DESeq2)
library(Biobase)
library(org.Hs.eg.db)
library(goseq)



counts_table <- read_csv("count_data/counts_table.csv")
select<-dplyr::select
d3<-counts_table%>%select(c(2,8:16))
colnames(d3)[2:10]<-c(paste("control",1:3,sep = "."), paste("treatment",1:6,sep = "."))

rownames(d3)<-counts_table$
d3<-d3%>%select(-c(1))

##1 
d4<-data.frame(HTSBasicFilter(d3, method ="sum",
                              cutoff.type = "value", cutoff = 10,
                              normalization = c("DESeq"))$filteredData)

##2
d5<-data.frame(HTSBasicFilter(d3, method ="rpkm",
                              cutoff.type = 2, cutoff = 0.5,
                              normalization = c("TMM"),length = counts_table$Length)$filteredData)



##3
d6<-data.frame(HTSBasicFilter(d3, method ="max",
                              cutoff.type = "quantile", cutoff = 0.3,
                              normalization = c("DESeq"),length = counts_table$Length
)$filteredData)



group<-factor(rep(c("A","B","C"),each =3))
c_d_1<-data.frame(sample.name = colnames(d4),group = group)



list_deseq_mat<-lapply(list(d4,d5,d6),function(x){DESeqDataSetFromMatrix(countData = x,colData = c_d_1,design = ~group)})


dds_list<-lapply(list_deseq_mat, function(x) { do.call("DESeq",args = list(x))})

y<-seq(0.01,0.1,by = 0.01)
vec_sum<-vec_max<-vec_rpkm<-numeric(10)
for(i in 1:length(y)){
  res_cutoff<-lapply(dds_list, function(x) { d1<- results(x); sum(d1$padj<y[i],na.rm = T) })
  vec_sum[i]<-res_cutoff[[1]]
  vec_max[i]<-res_cutoff[[2]]
  vec_rpkm[i]<-res_cutoff[[3]]
}
plot(y,vec_rpkm,typ = "b",col = 3 )
title(main = "green-rpkm,red-max,blue-sum ")
lines(y,vec_max,ty = "b",col = 2)
lines(y, vec_sum,typ = "b",col = 4)
abline(v= 0.05,col = "orange")

plotPCA(rlog(dds_list[[3]],blind = T), intgroup = c("group"))


mod_mat<-model.matrix(design(dds_list[[3]]), colData(dds_list[[3]]))
fac_lev<-lapply(LETTERS[1:3], function(x){ colMeans(mod_mat[dds_list[[3]]$group == x, ])})
si_rna_treated<- colMeans(mod_mat[dds_list[[3]]$group %in% c("B","C"),])
control_treated<-fac_lev[[1]]
res_df<-results(dds_list[[3]],contrast = si_rna_treated - control_treated,tidy = T)



res_df<-res_df %>% filter(padj<=0.05)


list_1<-list(up_regulated = res_df%>%filter(log2FoldChange>1) , down_regulated = res_df%>%filter(log2FoldChange<c(-1)))


# matching ensemble id to gene name


mart1<-useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

res1<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id',
            values = rownames(d3),mart = mart1)


# removing duplicated gene names
gene_list_with_name<- lapply(list_1, function(x){x1<-left_join(x,res1,by = c("row" = "ensembl_gene_id"));
                                                  x1<-x1 %>% filter(!duplicated(hgnc_symbol)); return(x1) } )


lapply(c(1,2), function(x) {
  sprintf("the number of %s genes found with p value = 0.05 is %i",names(gene_list_with_name)[x] ,nrow(gene_list_with_name[[x]]))
})

write.csv(gene_list_with_name[[1]], file = "upregulated.csv")
write.csv(gene_list_with_name[[2]], file = "downregulated.csv")

genes<-as.integer(results(dds_list[[3]])$padj<0.05)
not_na<-!is.na(genes)

names(genes)<-rownames(d6)
genes<-genes[not_na]

pwf<-nullp(genes,"hg19","ensGene")

go.mf<-goseq(pwf,"hg19","ensGene",test.cats = c("GO:MF"))
 
go.bp<-goseq(pwf,"hg19","ensGene",test.cats = c("GO:BP"))

write.csv(go.mf[p.adjust(go.mf$over_represented_pvalue,method = "BH")<0.05,],"over_represented_molecular function_go.csv")
write.csv(go.bp[p.adjust(go.bp$over_represented_pvalue,method = "BH")<0.05,],"over_represented_biological process_go.csv")












