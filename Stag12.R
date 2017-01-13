---
title: "TCGA clinical information"
author: "Matteo Cereda"
date: "01 MArch 2016"
output: html_document
---
  

fisher.test(matrix(c(3,129,110,(7827-242)), nr=2))

phyper(3, 129+110, 7827 - (110+129), 132, lower.tail=FALSE)

# config
source("config.R")

# SELECTION OF LOF ALTERATIONS ========

load("~/Lavoro/TCGA/total_table_muts_cnv.Rdata")
load("/Volumes/FC/DB/TCGA/01_03_2015/total_table_muts_cnv.Rdata")

ct = table(unique(total_table_muts_cnv[,1:2])$Cancer_type)

lof1 = as.data.frame(get_samples_driver_alterations(total_table_muts_cnv, "10274"), stringsAsFactors=F)
lof2 = as.data.frame(get_samples_driver_alterations(total_table_muts_cnv, "10735"), stringsAsFactors=F)

save(lof1,lof2, lof,file="Rdata/STAG1_STAG2.lof.Rdata")

lof = rbind(
  cbind(lof1, gene="STAG1"),
  cbind(lof2, gene="STAG2")
)

write.xlsx( lof, file="Results/STAG_1_2_LOF.xlsx","Samples", row.names = F,  showNA = F)

# table ========

s3=ddply(lof, .(Cancer_type), summarise, 
         n = length(unique(Sample)), n.f = length(unique(Sample[ gender!='MALE' & !is.na(gender)])), n.m = length(unique(Sample[ gender=='MALE' & !is.na(gender)])),
         STAG1 = length(unique(Sample[gene=="STAG1"])), STAG1_F = length(unique(Sample[gene=="STAG1" & gender!='MALE' & !is.na(gender)])), STAG1_M = length(unique(Sample[gene=="STAG1"  & gender=='MALE' & !is.na(gender)])),
         STAG2 = length(unique(Sample[gene=="STAG2"])), STAG2_F = length(unique(Sample[gene=="STAG2" &  gender!='MALE' & !is.na(gender)])), STAG2_M = length(unique(Sample[gene=="STAG2" & gender=='MALE' & !is.na(gender)]))
         
)

write.xlsx(s3, file="Results/Summary_sample_per_cancer_types.xlsx", row.names=F)

# BARPLOT ========

load("Rdata/STAG1_STAG2.lof.Rdata")

  # remove sample with mutations < 10%
# lof = lof[-which(lof$Sample=="TCGA-B0-4700-01"),]
# lof1 = lof1[-which(lof1$Sample=="TCGA-B0-4700-01"),]
# lof2 = lof2[-which(lof2$Sample=="TCGA-B0-4700-01"),]
# 
# lof$gender[which(lof$Sample=="TCGA-2G-AALR-01")]="MALE"
# lof1$gender[which(lof1$Sample=="TCGA-2G-AALR-01")]="MALE"
# lof2$gender[which(lof2$Sample=="TCGA-2G-AALR-01")]="MALE"
# 
# lof$simple_type="Truncating"
# lof$simple_type[which( lof$type=="Damaging mut." )] ="Damaging mut"
# 
# lof1$simple_type="Truncating"
# lof1$simple_type[which( lof1$type=="Damaging mut." )] ="Damaging mut"
# 
# lof2$simple_type="Truncating"
# lof2$simple_type[which( lof2$type=="Damaging mut." )] ="Damaging mut"

save(lof, lof1, lof2, file="Rdata/STAG1_STAG2.lof.Rdata")

load("Rdata/STAG1_STAG2.lof.Rdata")

x = table(lof1$Cancer_type, lof1$simple_type)
y = table(lof2$Cancer_type, lof2$simple_type)

x = melt(x); colnames(x)[1:2] = c("Var1","Var2"); x$p = x$value/ct[as.character(x$Var1)]*100
y = melt(y); colnames(y)[1:2] = c("Var1","Var2"); y$p = y$value/ct[as.character(y$Var1)]*100

# x$Var2 = factor(x$Var2, levels = c("Truncating mut.", "Damaging mut.", "More than 1 T/D muts", "Homozygous del.", "Heterozygous del. + trunc/dam mutation"), ordered = T)
# y$Var2 = factor(y$Var2, levels = c("Truncating mut.", "Damaging mut.", "More than 1 T/D muts", "Homozygous del.", "Heterozygous del. + trunc/dam mutation"), ordered = T)

x$Var2 = factor(x$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)
y$Var2 = factor(y$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)

x$Var1 = as.character(x$Var1)
y$Var1 = as.character(y$Var1)

tmp = ddply(x, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
x = rbind(x, cbind(Var1='ALL', tmp) )
x$value=as.numeric(x$value)
x$p = as.numeric(x$p)

tmp = ddply(y, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
y = rbind(y, cbind(Var1='ALL', tmp) )
y$value=as.numeric(y$value)
y$p = as.numeric(y$p)

x = subset(x, Var1%in%(y$Var1))

z = rbind(
  cbind(x, gene="STAG1"),
  cbind(y, gene="STAG2")
)
write.xlsx( z, file="Results/STAG_1_2_LOF.xlsx","LOF", row.names = F,  showNA = F, append = T)

# BARPLOT STAG2 M/F ========

m = with( subset(lof2, gender=="MALE"), table(Cancer_type, simple_type) )
f = with( subset(lof2, gender=="FEMALE"), table(Cancer_type, simple_type) )
m = melt(m); colnames(m)[1:2] = c("Var1","Var2"); m$p = m$value/ct[as.character(m$Var1)]*100
f = melt(f); colnames(f)[1:2] = c("Var1","Var2"); f$p = f$value/ct[as.character(f$Var1)]*100
m$Var2 = factor(m$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)
f$Var2 = factor(f$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)

m$Var1 = as.character(m$Var1)
f$Var1 = as.character(f$Var1)

tmp = ddply(m, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
m = rbind(m, cbind(Var1='ALL', tmp) )
m$value=as.numeric(m$value)
m$p = as.numeric(m$p)

tmp = ddply(f, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
f = rbind(f, cbind(Var1='ALL', tmp) )
f$value=as.numeric(f$value)
f$p = as.numeric(f$p)

z = rbind(z,
  cbind(m, gene="STAG2 (M)"),
  cbind(f, gene="STAG2 (F)")
)

m = with( subset(lof1, gender=="MALE"), table(Cancer_type, simple_type) )
f = with( subset(lof1, gender=="FEMALE"), table(Cancer_type, simple_type) )
m = melt(m); colnames(m)[1:2] = c("Var1","Var2"); m$p = m$value/ct[as.character(m$Var1)]*100
f = melt(f); colnames(f)[1:2] = c("Var1","Var2"); f$p = f$value/ct[as.character(f$Var1)]*100
m$Var2 = factor(m$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)
f$Var2 = factor(f$Var2, levels = c("Truncating", "Damaging mut"), ordered = T)
m$Var1 = as.character(m$Var1)
f$Var1 = as.character(f$Var1)

tmp = ddply(m, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
m = rbind(m, cbind(Var1='ALL', tmp) )
m$value=as.numeric(m$value)
m$p = as.numeric(m$p)

tmp = ddply(f, .(Var2), summarize, value=sum(value)); tmp$p = (tmp$value/7828)*100
f = rbind(f, cbind(Var1='ALL', tmp) )
f$value=as.numeric(f$value)
f$p = as.numeric(f$p)

m = subset(m, Var1%in%unique(z$Var1))
f = subset(f, Var1%in%unique(z$Var1))


z = rbind(z,
          cbind(m, gene="STAG1 (M)"),
          cbind(f, gene="STAG1 (F)")
)

y = table(lof2$Cancer_type, lof2$simple_type)
y = melt(y); colnames(y)[1:2] = c("Var1","Var2"); y$p = y$value/ct[as.character(y$Var1)]*100

tmp = ddply( subset(y, Var1!="ALL"), .(Var1), summarise, n=sum(p))
z$Var1 = factor( z$Var1, levels= c(as.character(tmp$Var1[order(tmp$n)]), "ALL"))
z$Var2 = factor(z$Var2, levels = c("Truncating","Damaging mut"), ordered = T)
z$gene = factor(z$gene, levels=c("STAG2 (M)", "STAG2 (F)",'STAG2',"STAG1 (M)", "STAG1 (F)",'STAG1'))
pdf(file="Results/STAG_1_2_LOF.pdf", h=18, w=12)
ggplot(arrange(z, Var2), aes(x=Var1,y=p, fill=Var2))+geom_bar(stat="identity")+theme_boss()+ylab('')+xlab('')+
  theme(axis.text.x = element_text(angle=90), strip.background=element_blank())+facet_grid(gene~.)+scale_y_continuous(limits=c(0,12), breaks=c(0,6,12))#+scale_fill_manual(values=c('blue','red'))+
dev.off()

write.xlsx( subset(z, gene%in%c("STAG2 (M)","STAG2 (F)")), file="Results/STAG_1_2_LOF.xlsx","LOF STAG2 M,F", row.names = F,  showNA = F, append = T)

x = read.xlsx("Results/STAG_1_2_LOF.xlsx", "One mutation")
x = subset(x, gene=="STAG2")
x$patient  = substr(x$sample, 1, 12)
y =x
load("Rdata/STAG1_STAG2.gender.Rdata")
y$gender  = x[match(y$patient, x[,1]),2]
write.xlsx( y, file="Results/STAG_1_2_LOF.xlsx","One mut. STAG2 M,F", row.names = F,  showNA = F, append = T)
x = read.xlsx("Results/STAG_1_2_LOF.xlsx", "One mutation")
x = subset(x, gene=="STAG1")
x$patient  = substr(x$sample, 1, 12)
y =x
load("Rdata/STAG1_STAG2.gender.Rdata")
y$gender  = x[match(y$patient, x[,1]),2]
y$MeanMutFreq = as.numeric(y$MeanMutFreq)
write.xlsx( y, file="Results/STAG_1_2_LOF.xlsx","One mut. STAG1 M,F", row.names = F,  showNA = F, append = T)

write.xlsx( subset(z, gene%in%c("STAG1 (M)", "STAG1 (F)")), file="Results/STAG_1_2_LOF.xlsx","One mut. STAG1 M,F", row.names = F,  showNA = F, append = T)

## Figure 1D =====

load("Rdata/STAG1_STAG2.lof.Rdata")
load("Rdata/barplot_order.Rdata")

m = with( subset(lof2, gender=="MALE"), table(Cancer_type, simple_type) ); sm = apply(m,1,sum)
f = with( subset(lof2, gender=="FEMALE"), table(Cancer_type, simple_type) ); sf = fapply(f,1,sum)
m = melt(m); colnames(m)[1:2] = c("Var1","Var2"); m$p = m$value/sm*100; m$label = paste0(m$Var1, " (", sm[m$Var1],")")
f = melt(f); colnames(f)[1:2] = c("Var1","Var2"); f$p = f$value/sf*100; f$label = paste0(f$Var1, " (", sf[f$Var1],")")
m$Var2 = as.character(m$Var2)
f$Var2 = as.character(f$Var2)
m$Var1 = as.character(m$Var1)
f$Var1 = as.character(f$Var1)

z = rbind(cbind(m, gene="STAG1 (M, 62)"),
          cbind(f, gene="STAG1 (F, 70)")
)
z = z[rev(order(z[,2])),]

z$Var1 = factor(z$Var1, levels=plot_order)
z$Var2 = factor(z$Var2, levels = c("Damaging mut","Truncating"), ordered = T)

write.xlsx( z, file="Results/STAG_2_LOF_male_female.xlsx","LOF STAG2 M,F", row.names = F,  showNA = F, append = T)


pdf(file="Results/STAG_2_LOF_male_female.pdf", h=6, w=6)
ggplot(z, aes(x=Var1, y = p, fill=Var2))+geom_bar(stat='identity')+facet_wrap(~gene, ncol=1)+theme_classic()+theme(axis.text.x = element_text(angle=90))
dev.off()

pdf(file="Results/STAG_2_LOF_male_female.pie.pdf", h=6, w=12)
par(mfrow=c(1,2))
m = with( subset(lof2, gender=="MALE"), table(Cancer_type, simple_type) ); sm = sum(m)
uno = apply(m,2, sum)/sm; names(uno) = as.character(round(uno,3))
pie(uno, clockwise = T)
f = with( subset(lof2, gender=="FEMALE"), table(Cancer_type, simple_type) ); sf = sum(f)
uno = apply(f,2, sum)/sf; names(uno) = as.character(round(uno,3))
pie(uno, clockwise = T)
dev.off()

# Coocurrence of mutations ========

colof1 = subset(lof1,  no_TRUNC_muts>1 | no_NTDam_muts>1 | (no_TRUNC_muts==1 & no_NTDam_muts==1) )
colof2 = subset(lof2,  no_TRUNC_muts>1 | no_NTDam_muts>1 | (no_TRUNC_muts==1 & no_NTDam_muts==1) )

x = rbind()
for(i in 1:nrow(colof1)){
  print(i)
  load(paste0("/Volumes/FC/DB/TCGA/01_03_2015/",colof1[i,1],"/Tumor/Somatic_Mutations/",colof1[i,1],"_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust_LOFGOF.Rdata"))
  tmp = subset(somatic_mutations[[ colof1$Sample[i] ]], !duplicated(key) & symbol=="STAG1" )
  x = rbind(x, cbind('cancer'=colof1[i,1], 'sample'=colof1$Sample[i], tmp[,c('Chr','Start','End', 'Ref', 'Alt',"ExonicFunc.refGene","MeanMutFreq",'Sequence_Source')] ) )
}
# 
x$key = paste0(x$sample,".",x$Chr,".",x$End,".",x$Alt,".",x$Ref, sep="", coll="")
x$duplicated = duplicated(x$key)

y = rbind()
for(i in 1:nrow(colof2)){
  print(i)
  load(paste0("/Volumes/FC/DB/TCGA/01_03_2015/",colof2[i,1],"/Tumor/Somatic_Mutations/",colof2[i,1],"_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust_LOFGOF.Rdata"))
  tmp = subset(somatic_mutations[[ colof2$Sample[i] ]], !duplicated(key) & symbol=="STAG2" )
  y = rbind(y, cbind('cancer'=colof2[i,1], 'sample'=colof2$Sample[i], tmp[,c('Chr','Start','End', 'Ref', 'Alt',"ExonicFunc.refGene","MeanMutFreq",'Sequence_Source')] ) )
}

z = rbind(
  cbind(x[,c('cancer','sample','Chr','End','Ref','Alt','ExonicFunc.refGene','MeanMutFreq')], gene="STAG1"),
  cbind(y[,c('cancer','sample','Chr','End','Ref','Alt','ExonicFunc.refGene','MeanMutFreq')], gene="STAG2")
)

z$type="Truncating mut."
z$type[which( z$ExonicFunc.refGene%in%c('splicing',"nonsynonymous SNV") ) ]="Damaging mut."
z$key = paste0(z$sample,"  ",z$Chr,".",z$End,".",z$Alt,".",z$Ref, sep="", coll="")

write.xlsx(z[,1:9], file="Results/Stag2_LOF_Stag1_LOF.cooccurrence.xlsx", "Muts",row.names = F)

# Frequenza di quelle con una sola mutazione ========

colof1 = subset(lof,  gene=="STAG1" & type%in%c("Damaging mut.","Truncating mut."))
colof2 = subset(lof,  gene=="STAG2" & type%in%c("Damaging mut.","Truncating mut."))

x = rbind()
cc = ""
for(i in 1:nrow(colof1)){
  print(i)
  if(colof1[i,1]!=cc){
  load(paste0("/Volumes/FC/DB/TCGA/01_03_2015/",colof1[i,1],"/Tumor/Somatic_Mutations/",colof1[i,1],"_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust_LOFGOF.Rdata"))
  }
  cc = colof1[i,1]
  tmp = subset(somatic_mutations[[ colof1$Sample[i] ]], !duplicated(key) & symbol=="STAG1" )
  x = rbind(x, cbind('cancer'=colof1[i,1], 'sample'=colof1$Sample[i], tmp[,c('Chr','Start','End', 'Ref', 'Alt',"ExonicFunc.refGene","MeanMutFreq",'Sequence_Source')] ) )
}

y = rbind()
for(i in 1:nrow(colof2)){
  print(i)
  if(colof2[i,1]!=cc){
    load(paste0("/Volumes/FC/DB/TCGA/01_03_2015/",colof2[i,1],"/Tumor/Somatic_Mutations/",colof2[i,1],"_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust_LOFGOF.Rdata"))
  }
  cc = colof2[i,1]
  tmp = subset(somatic_mutations[[ colof2$Sample[i] ]], !duplicated(key) & symbol=="STAG2" )
  y = rbind(y, cbind('cancer'=colof2[i,1], 'sample'=colof2$Sample[i], tmp[,c('Chr','Start','End', 'Ref', 'Alt',"ExonicFunc.refGene","MeanMutFreq",'Sequence_Source')] ) )
}

z = rbind(
  cbind(gene="STAG1", x[,c('cancer','sample','Chr','End','Ref','Alt','ExonicFunc.refGene','MeanMutFreq')]),
  cbind(gene="STAG2", y[,c('cancer','sample','Chr','End','Ref','Alt','ExonicFunc.refGene','MeanMutFreq')])
)

z$type="Truncating mut."
z$type[which( z$ExonicFunc.refGene%in%c('splicing',"nonsynonymous SNV") ) ]="Damaging mut."

write.xlsx( z, file="Results/STAG_1_2_LOF.xlsx","One T/D mut.", row.names = F,  showNA = F, append = T)
# ggplot(z, aes(x=MeanMutFreq))+geom_histogram(binwidth = 0.05,fill="grey80",col="black")+theme_boss()+facet_grid(~gene)+coord_cartesian()+scale_x_continuous(limits = c(0,1))

# EXPRESSION STAG1 STAG2 ========

# df = total_table_muts_cnv
# lof1 = as.data.frame(get_samples_driver_alterations(df, "10274"), stringsAsFactors=F)
# lof2 = as.data.frame(get_samples_driver_alterations(df, "10735"), stringsAsFactors=F)
# 
# 
# rt1 = subset(total_rna, Entrez=="10274")
# rt2 = subset(total_rna, Entrez=="10735")
# 
# rt1$stag1_m="WT"
# rt1$stag2_m="WT"
# rt1$stag1_m[which(rt1$Sample%in%lof1[,2])]="MUT"
# rt1$stag2_m[which(rt1$Sample%in%lof2[,2])]="MUT"
# 
# rt2$stag1_m="WT"
# rt2$stag2_m="WT"
# rt2$stag1_m[which(rt2$Sample%in%lof1[,2])]="MUT"
# rt2$stag2_m[which(rt2$Sample%in%lof2[,2])]="MUT"
# 
# save(rt1, file="~/OAC/rna_stag1.Rdata")
# save(rt2, file="~/OAC/rna_stag2.Rdata")

load("Rdata/rna_stag1.Rdata")
load("Rdata/rna_stag2.Rdata")
load("Rdata/STAG1_STAG2.lof.Rdata")

rt2$type = lof2[match(rt2$Sample, lof2$Sample),"type"]
rt2$type[which(is.na(rt2$type))] = "WT"

x = subset(rt2, stag2_m=="MUT")
x$type = "Mutated"

x = rbind(x, rt2)

x$type = factor(x$type, levels=c('Mutated',"Damaging mut.","Truncating mut.","More than 1 T/D muts","Homozygous del.","Heterozygous del. + trunc/dam mutation","WT"))
t = table(x$type)
x$labels = paste0("WT (", t["WT"], ")") 
x$labels[ which(x$type=="Mutated")] = paste0("Mutated (", t["Mutated"], ")")
x$labels[ which(x$type%in%c("Truncating mut.", 
                            "Homozygous del.", 
                            "Heterozygous del. + trunc/dam mutation",
                            "More than 1 T/D muts")) ] = paste0("Truncating (", sum(t[c("Truncating mut.", 
                                                                                                         "Homozygous del.", 
                                                                                                         "Heterozygous del. + trunc/dam mutation",
                                                                                        "More than 1 T/D muts")]), ")")
x$labels[ which(x$type=="Damaging mut.")] = paste0("Damaging mut. (", t["Damaging mut."], ")")
x$labels = factor(x$labels, levels=c("Mutated (133)", "Damaging mut. (71)", "Truncating (62)","WT (7695)"))

pdf(file="Results/STAG2_TPM.pdf",h=8,w=6)
ggplot(x, aes(y=TPM,x=labels))+geom_boxplot()+theme_boss()+xlab("")+
  ggtitle(paste0("P-value =",with(rt2, round(wilcox.test(TPM[stag2_m=="MUT"], TPM[stag2_m=="WT"])$p.value,10)),
                 "\nP-value (dam) =",with(rt2, round(wilcox.test(TPM[type=="Damaging mut."], TPM[type=="WT"])$p.value,3)),
                 "\nP-value (tr) =",with(x, round(wilcox.test(TPM[labels=="Truncating (62)"], TPM[labels=="WT (7695)"])$p.value,21)),
                 "\nP-value (>1) =",with(rt2, round(wilcox.test(TPM[type=="More than 1 T/D muts"], TPM[type=="WT"])$p.value,3)),
                 "\nTruncating = trunc muts + Homo del + \n( Het del + trunc/dam muts) +(More than 1 T/D muts)"
                 )
          )+theme(axis.text.x=element_text(angle=45, hjust = 1))
dev.off()
write.xlsx(rt2, file="Results/STAG2_TPM.xlsx", row.names = F)

# STAG1

rt1$type = lof2[match(rt1$Sample, lof2$Sample),"type"]
rt1$type[which(is.na(rt1$type))] = "WT"

x = subset(rt1, stag2_m=="MUT")
x$type = "Mutated"

x = rbind(x, rt1)

x$type = factor(x$type, levels=c('Mutated',"Damaging mut.","Truncating mut.","More than 1 T/D muts","Homozygous del.","Heterozygous del. + trunc/dam mutation","WT"))
t = table(x$type)
x$labels = paste0("WT (", t["WT"], ")") 
x$labels[ which(x$type=="Mutated")] = paste0("Mutated (", t["Mutated"], ")")
x$labels[ which(x$type%in%c("Truncating mut.", 
                            "Homozygous del.", 
                            "Heterozygous del. + trunc/dam mutation",
                            "More than 1 T/D muts")) ] = paste0("Truncating (", sum(t[c("Truncating mut.", 
                                                                                        "Homozygous del.", 
                                                                                        "Heterozygous del. + trunc/dam mutation",
                                                                                        "More than 1 T/D muts")]), ")")
x$labels[ which(x$type=="Damaging mut.")] = paste0("Damaging mut. (", t["Damaging mut."], ")")
x$labels = factor(x$labels, levels=c("Mutated (133)", "Damaging mut. (71)", "Truncating (62)","WT (7695)"))

pdf(file="Results/STAG1_when_STAG2_mut_TPM.pdf",h=8,w=6)
ggplot(x, aes(y=TPM,x=labels))+geom_boxplot()+theme_boss()+xlab("")+
  ggtitle(paste0("P-value =",with(rt1, round(wilcox.test(TPM[stag2_m=="MUT"], TPM[stag2_m=="WT"])$p.value,3)),
                 "\nP-value (dam) =",with(rt1, round(wilcox.test(TPM[type=="Damaging mut."], TPM[type=="WT"])$p.value,3)),
                 "\nP-value (tr) =",with(x, round(wilcox.test(TPM[labels=="Truncating (62)"], TPM[labels=="WT (7695)"])$p.value,3)),
                 "\nP-value (>1) =",with(rt1, round(wilcox.test(TPM[type=="More than 1 T/D muts"], TPM[type=="WT"])$p.value,3)),
                 "\nTruncating = trunc muts + Homo del + \n( Het del + trunc/dam muts) +(More than 1 T/D muts)"
                 
  )
  )+theme(axis.text.x=element_text(angle=45, hjust = 1))
dev.off()
write.xlsx(rt1, file="Results/STAG1_TPM.xlsx", row.names = F)







p1=ggplot(rt1, aes(y=TPM,x=stag1_m))+geom_boxplot()+theme_boss()+ggtitle(paster("TPM STAG1 - wilcoxon p =",
                                                                             with(rt1, round(wilcox.test(TPM[stag1_m=="MUT"], TPM[stag1_m=="WT"])$p.value,3))))+xlab("STAG1")

with(rt2, wilcox.test(TPM[stag2_m=="MUT"], TPM[stag2_m=="WT"])$p.value)
p2=
  ggplot(rt2, aes(y=TPM,x=stag2_m))+geom_boxplot()+theme_boss()+ggtitle("TPM STAG2 - wilcoxon p = 2.9x10-8")+xlab("STAG2")

pdf(file="Results/Expression_STAG1_STAG2.pdf", h=6, w=8)
grid.arrange(p1,p2, nrow=1)
dev.off()

stats= data.frame( 
  "STAG1_MUT"=as.matrix(summary( subset(rt1, stag1_m=="MUT")$TPM )),
  "STAG1_WT"=as.matrix(summary( subset(rt1, stag1_m=="WT")$TPM )),
  "STAG2_MUT"=as.matrix(summary( subset(rt2, stag2_m=="MUT")$TPM )),
  "STAG2_WT"=as.matrix(summary( subset(rt2, stag2_m=="WT")$TPM ))
)
write.xlsx(stats, file="Results/Expression_STAG1_STAG2.xlsx","stats")
write.xlsx(rt1, file="Results/Expression_STAG1_STAG2.xlsx", "STAG1", append=T)
write.xlsx(rt2, file="Results/Expression_STAG1_STAG2.xlsx", "STAG2", append=T)

p1=ggplot(data=r, aes(x=lab, y=TPM))+geom_boxplot()+ theme_boss()+xlab("STAG2")+ggtitle(paster("TPM STAG1 - wilcoxon p =",
round(wilcox.test(r$TPM[r$STAG2=="MUT"], r$TPM[r$STAG2=="WT"])$p.value,3)))

p2=ggplot(data=r, aes(x=lab, y=TPM))+geom_boxplot()+ theme_boss()+xlab("STAG1")+ggtitle(paster("TPM STAG2 - wilcoxon p =",
round(wilcox.test(r$TPM[r$STAG1=="MUT"], r$TPM[r$STAG1=="WT"])$p.value,3)))

pdf(file="Results/TPM_stag1_stag2.pdf", h=8, w=8)
grid.arrange(p1,p2,nrow=1)
dev.off()

# CO-MUTATED STAG1 STAG2 ========

df = subset(total_table_muts_cnv , !is.na(no_TRUNC_muts) | !is.na(no_NTDam_muts) | (CNV_type=="Loss" & Copy_number==0) )

lof1 =subset(df, Entrez=="10274", select = c("Cancer_type","Sample"))
lof2 =subset(df, Entrez=="10735", select = c("Cancer_type","Sample"))

pdf(file="Results/venn_stag1_stag2.pdf",w=4,h=4)
venn(list(stag1=lof1[,2],stag2=lof2[,2]))
dev.off()

u = data.frame(cancer_type=NA, sample=union(lof1[,2],lof2[,2]), stag1_mut=F,stag2_mut=F)
u$stag1_mut[which(u$sample%in%lof1[,2])]=T
u$stag2_mut[which(u$sample%in%lof2[,2])]=T
u$cancer_type[which(u$sample%in%lof1[,2])]=lof1[,1]
u$cancer_type[which(u$sample%in%lof2[,2])]=lof2[,1]

write.table(u, file="Results/stag1_stag2_samples.tsv", row.names = F)

x = subset(df, Cancer_type %in% c("BLCA","KIRC", "SKCM", "UCEC") )
y = ddply(x, .(Cancer_type,Sample), summarize, no_NSI_muts = sum(no_NSI_muts, na.rm=T), cnv=sum(!is.na(CNV_type)), f.cnv=sum(!is.na(CNV_type))/18956, .progress="text")
y$stag1_2_mut=F
y$stag1_2_mut[which(y$Sample %in% subset(u, stag1_mut & stag2_mut)$sample)]=T

y=y[order(y[,1],y[,3]),]

p1 = ggplot(y, aes(x=Cancer_type, y=no_NSI_muts))+geom_boxplot()+ggtitle("Number of nonsilent mutations")+geom_point(data=subset(y,stag1_2_mut), aes(x=Cancer_type, y=no_NSI_muts,fill=stag1_2_mut),size =4, shape=21)+ theme_boss()
p2 = ggplot(y, aes(x=Cancer_type, y=cnv))+geom_boxplot()+ggtitle("Number of CNV genes")+geom_point(data=subset(y,stag1_2_mut), aes(x=Cancer_type, y=cnv,fill=stag1_2_mut),size =4, shape=21)+ theme_boss()
pdf(file="~/Lavoro/TCGA/check_samples_stag1_2_mt.pdf",h=8,w=8)
grid.arrange(p1,p2,nrow=2)
dev.off()
ddply(y, .(Cancer_type), function(x) s(x$no_NSI_muts))
ddply(y, .(Cancer_type), function(x) s(x$cnv))

# Exploratory comorbity analysis
fisher.test(matrix(c( 4, 110, 129, 7561 ),nr=2,nc=2,byrow=T))

p.STAG1 = 114/7561
p.STAG2 = 133/7561
expected = p.STAG1*p.STAG2*7561

# fisher.test(matrix(c( 1, 110, 186, 365 ),nr=2,nc=2,byrow=T))


# Fisher's Exact Test for Count Data
# 
# data:  x
# p-value = 0.1297
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.5620566 5.7464081
# sample estimates:
# odds ratio 
# 2.131069 

x = subset(df, Sample%in%subset(u, stag1_mut & stag2_mut)$sample & Entrez%in%c("10274","10735"))

load("/Volumes/FC/DB/TCGA/01_03_2015/BLCA/Tumor/Somatic_Mutations/BLCA_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
y = cbind(cancer="BLCA",subset(somatic_mutations[["TCGA-R3-A69X-01"]], entrez_19014 %in% c("10274","10735") ))

load("/Volumes/FC/DB/TCGA/01_03_2015/KIRC/Tumor/Somatic_Mutations/KIRC_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
y = rbind(y, cbind(cancer="KIRC",subset(somatic_mutations[["TCGA-B0-4700-01"]], entrez_19014 %in% c("10274","10735") )))

load("/Volumes/FC/DB/TCGA/01_03_2015/SKCM/Tumor/Somatic_Mutations/SKCM_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
y = rbind(y, cbind(cancer="SKCM",subset(somatic_mutations[["TCGA-EE-A2A2-06"]], entrez_19014 %in% c("10274","10735") )))

load("/Volumes/FC/DB/TCGA/01_03_2015/UCEC/Tumor/Somatic_Mutations/UCEC_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF_Annovar_dbNSFP_gAnn_SampleFilters_NonSilent_MutationFilters_OncodriveClust.Rdata")
y = rbind(y, cbind(cancer="UCEC",subset(somatic_mutations[["TCGA-A5-A0GP-01"]], entrez_19014 %in% c("10274","10735") )))

load("/Volumes/FC/DB/TCGA/01_03_2015/KIRC/Tumor/Somatic_Mutations/KIRC_concatenated_deduplicated_MutExpCNV_1SamplePerPatient_MMF.Rdata")

write.xlsx(y, file="Results/comutated_samples.STAG1_STAG2.mutations.xlsx", showNA = F, row.names = F )
# co-mutated ALL ENTREZ ====

df = subset(total_table_muts_cnv , !is.na(no_TRUNC_muts) | !is.na(no_NTDam_muts) | (CNV_type=="Loss" & Copy_number==0) )
rm(total_table_muts_cnv)
df = as.data.frame(df, stringsAsFactors=F)

entrez = as.list(unique(df$Entrez))
entrez = lapply_pb(entrez, function(x,d) subset(df, Entrez==x, select = c("Sample"))[,1], d=df)
names(entrez) = unique(df$Entrez)

get_comutated =  function(x,y) length(intersect(x,y))
e = unique(df$Entrez)

int = vector("list", length(entrez))
for(i in 1:length(entrez)){
  print(i)
  int[[i]] = sapply(entrez, get_comutated, y= subset(df, Entrez==e[i], select = c("Sample"))[,1])
}
names(int) = unique(df$Entrez)

a = do.call("rbind", int)
b=a
b[upper.tri(b)]=NA
comut = b

colnames(comut) = as.character(unique(df$Entrez))
rownames(comut) = as.character(unique(df$Entrez))

save(comut, file="Rdata/comut_matrix_driver_alterations.Rdata")
load("Rdata/comut_matrix_driver_alterations.Rdata")
load("Rdata/tmp_entrez.Rdata")

comut[ row(comut)==col(comut) ] = NA

a = na.omit(c(comut)); quantile(a, seq(0,1,.05))

heatmap(comut)


stag1 = int[["10274"]]
names(stag1)= unique(df$Entrez)
stag1 = stag1[which(names(stag1)!="10274")]
s1 = data.frame(id="STAG1", genes=names(stag1), value=stag1)
pem = sum(stag1>4)/(length(stag1))
p1=ggplot(s1, aes(x=id, y=value))+geom_boxplot()+theme_boss()+xlab("")+geom_hline(yintercept = 4, col="red")+ylab('Number of samples')+
  ggtitle(paste0("Number of samples where STAG1\nis comutated with another gene\np-value (STAG1,STAG2) = ",round(pem,3)))
s(stag1)

stag2 = int[["10735"]]
names(stag2)= unique(df$Entrez)
stag2 = stag2[which(names(stag2)!="10735")]
s2 = data.frame(id="STAG2", genes=names(stag2), value=stag2)
pem=sum(stag2>4)/(length(stag2))
p2=ggplot(s2, aes(x=id, y=value))+geom_boxplot()+theme_boss()+xlab("")+geom_hline(yintercept = 4, col="red")+ylab('Number of samples')+
  ggtitle(paste0("Number of samples where STAG2 is\ncomutated with another gene\np-value (STAG2,STAG1) = ",round(pem,3)))
s(stag2)

pdf(file="Results/comutated_samples.STAG1_STAG2.pdf",h=6, w=8)
grid.arrange(p1,p2, ncol=2)
dev.off()

x = rbind(s1,s2)
write.xlsx(x, file="Results/comutated_samples.STAG1_STAG2.xlsx", "comutated samples",row.names=F)
write.xlsx(as.matrix(s(stag1)), file="Results/comutated_samples.STAG1_STAG2.xlsx", "STAG1 dist.",appen=T, row.names=F)
write.xlsx(as.matrix(s(stag2)), file="Results/comutated_samples.STAG1_STAG2.xlsx", "STAG2 dist.",appen=T, row.names=F)

