#environment
library(tidyverse)
library(tidytext)
library(vegan)
library(data.table)
library(ROCR)
library(lmerTest)
library(pheatmap)
library(PResiduals)
library(grid) 
library(RColorBrewer)

#data transformation at genus level
#OTUs a list of OTU table that is defined in file Diversity.R and saveRDs
OTUs<-readRDS("OTUs")
metadata<-read.csv("metadata.csv",header=T)
Genus.CZM<-zCompositions::cmultRepl(t(OTUs$TaxaSummary$Genus),
                                   method="CZM",
                                   output="p-counts",
                                   suppress.print=TRUE) %>% t() #0-replaced
Genus.CLR<-apply(log2(Genus.CZM), 2, function(x)x-mean(x)) #CLR-normalized

data<-t(Genus.CLR) %>% as.data.frame() %>% rownames_to_column("id")
data$age <- metadata$age[match(data$id, metadata$id)]
data$study<- metadata$study[match(data$id, metadata$id)]
data<-data %>% as.data.frame() %>% select(age, study,everything()) %>% 
  mutate(age=factor(age, levels=c("young","old"))) 

#liner mixed effect analysis
genus.lmer.results<-data.frame() 
for (i in 4:length(data)){
  taxa<-colnames(data)[i]
  mod <- lmerTest:::lmer(data[,i] ~ data$age+(1|study),data=data)
  meandifference<- summary(mod)$coefficients[2,1]
  pval <- summary(mod)$coefficients[2,5]
  cf<-confint(mod)
  cf_low <- cf[4,1]
  cf_high <- cf[4,2]
  genus.lmer.results <- rbind(genus.lmer.results, data.frame(Taxa=taxa, 
                                                      meanofdifference=meandifference,
                                                      CI_low=cf_low,
                                                      CI_high=cf_high,
                                                      Pval = pval))
}
genus.lmer.results$FDR_BH <- p.adjust(genus.lmer.results$Pval, method = 'BH')

P.genus.filter<-genus.lmer.results %>% 
  mutate(Significant=case_when(
    FDR_BH<0.05 & meanofdifference >0 ~ "increased in old",
    FDR_BH<0.05 & meanofdifference<0 ~ "increased in young",
    TRUE~"not significant"
  ))%>% filter(FDR_BH<0.05) %>% separate(Taxa, 
         c("Kingdom",
           "Phylum",
           "Class",
           "Order",
           "Family",
           "Genus"
         ), 
         sep=";", 
         remove=FALSE)

#rename lmer results description content
P.genus.filter$Des.taxa<-paste(P.genus.filter$Phylum,P.genus.filter$Genus,sep=";")
P.genus.filter[P.genus.filter=="p__Bacteroidetes;g__"]<-"p__Bacteroidetes;g__unclassified"
P.genus.filter[P.genus.filter=="p__Proteobacteria;g__"]<-"p__Proteobacteria;g__unclassified"
P.genus.filter[P.genus.filter=="p__TM7;g__"]<-"p__TM7;g__unclassified"

#visualization by ggplot2
level<-P.genus.filter %>% arrange(desc(Significant)) %>% pull(Des.taxa)
P.genus<-P.genus.filter %>% mutate(Des.taxa=factor(Des.taxa, levels=level)) %>% ungroup() %>%
  ggplot(aes(x=meanofdifference, y=Des.taxa, color=Significant)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
  geom_point() +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1),axis.text = element_text(color="black"))+
  scale_color_manual(values=c("firebrick3","dodgerblue3","black"))+
  xlab("log2FC (95% CI)")+
  ylab("Taxa at genus level")+
  ggtitle("Taxa enriched by age")

tiff('Fig_final/taxa.enrich.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')
P.genus
dev.off()

#data transformation at phylum level
Phylum.CZM<-zCompositions::cmultRepl(t(OTUs$TaxaSummary$Phylum),
                                     method="CZM",
                                     output="p-counts",
                                     suppress.print=TRUE) %>% t() #0-replaced
Phylum.CLR<-apply(log2(Phylum.CZM), 2, function(x)x-mean(x)) #CLR-normalized

data<-t(Phylum.CLR) %>% as.data.frame() %>% rownames_to_column("id")
data$age <- metadata$age[match(data$id, metadata$id)]
data$study<- metadata$study[match(data$id, metadata$id)]
data<-data %>% as.data.frame() %>% select(age, study,everything()) %>% 
  mutate(age=factor(age, levels=c("young","old"))) 

phylum.lmer.results<-data.frame() 
for (i in 4:length(data)){
  taxa<-colnames(data)[i]
  mod <- lmerTest:::lmer(data[,i] ~ data$age+(1|study),data=data)
  meandifference<- summary(mod)$coefficients[2,1]
  pval <- summary(mod)$coefficients[2,5]
  cf<-confint(mod)
  cf_low <- cf[4,1]
  cf_high <- cf[4,2]
  phylum.lmer.results <- rbind(phylum.lmer.results, data.frame(Taxa=taxa, 
                                                               meanofdifference=meandifference,
                                                               CI_low=cf_low,
                                                               CI_high=cf_high,
                                                               Pval = pval))
}
phylum.lmer.results$FDR_BH <- p.adjust(phylum.lmer.results$Pval, method = 'BH')

#export results
write.csv(phylum.lmer.results,"phylum.lmer.results.csv")
write.csv(genus.lmer.results,"genus.lmer.results.csv")

#data transformation for pathway dataset and lmer analysis for pathway
data.pw<-PW.RA %>% column_to_rownames(.,var="pathway") %>% t(.) %>% 
  as.data.frame() %>% rownames_to_column("id")

data.pw<-zCompositions::cmultRepl(t(PW.filter),
                                    method="CZM",
                                    output="p-counts",
                                    suppress.print=TRUE) %>% t() #0-replaced
data.pw<-apply(log2(data.pw), 2, function(x)x-mean(x)) #CLR-normalized
data.pw<-t(data.pw)%>% as.data.frame() %>% rownames_to_column("id")
data.pw$age <- metadata$age[match(data.pw$id, metadata$id)]
data.pw$study<- metadata$study[match(data.pw$id, metadata$id)]
data.pw<-data.pw %>% as.data.frame() %>% select(age, study,everything()) %>% 
  mutate(age=factor(age, levels=c("young","old"))) 

pw.lmer.results<-data.frame() 
for (i in 4:length(data.pw)){
  pw<-colnames(data.pw)[i]
  mod <- lmerTest:::lmer(data.pw[,i] ~ data.pw$age+(1|study),data=data.pw)
  meandifference<- summary(mod)$coefficients[2,1]
  pval <- summary(mod)$coefficients[2,5]
  cf<-confint(mod)
  cf_low <- cf[4,1]
  cf_high <- cf[4,2]
  pw.lmer.results <- rbind(pw.lmer.results, data.frame(Pathway=pw, 
                                                             meanofdifference=meandifference,
                                                             CI_low=cf_low,
                                                             CI_high=cf_high,
                                                             Pval = pval))
}

#adjust p value by BH method
pw.lmer.results$FDR_BH <- p.adjust(pw.lmer.results$Pval, method = 'BH')

pw.lmer.results$des<-PW.des$description[match(pw.lmer.results$Pathway, PW.des$feature)]

P.pw.filter<-pw.lmer.results %>% 
  mutate(Significant=case_when(
    FDR_BH<0.05 & meanofdifference >0 ~ "increased in old",
    FDR_BH<0.05 & meanofdifference<0 ~ "increased in young",
    TRUE~"not significant"
  ))%>% filter(FDR_BH<0.05)

level<-P.pw.filter %>% arrange(desc(Significant)) %>% pull(des)

P.pathway<-P.pw.filter %>% mutate(Des.pw=factor(des, levels=level)) %>% ungroup() %>%
  ggplot(aes(x=meanofdifference, y=Des.pw, color=Significant)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
  geom_point() +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1),axis.text = element_text(color="black"))+
  scale_color_manual(values=c("firebrick3","dodgerblue3","black"))+
  xlab("log2FC (95% CI)")+
  ylab("Pathway")+
  ggtitle("Pathway enriched by age")

tiff('Fig_final/pw.enrich.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')
P.pathway
dev.off()

#spearman cor analysis genus vs gpr41/43
#import dataset of relative abundance of taxa and gpr41/43 gene expression
cor.genus<-read.csv("taxa_vs_gprs.csv",header=T,row.names = 1)

#genus vs gpr41
cor.genus.result.gpr41<-list()
for (i in 3:length(cor.genus)){
  genus<-colnames(cor.genus)[i]
  mod <- cor.test(cor.genus[,i], cor.genus$Gpr41, method=c("spearman"))
  pval <- mod$p.value
  r<-mod$estimate
  cor.genus.result.gpr41 <- rbind(cor.genus.result.gpr41, data.frame(genus=genus, 
                                                       pval=pval,
                                                       r=r))
}

#adjust for p value by BH method
cor.genus.result.gpr41$FDR_BH <- p.adjust(cor.genus.result.gpr41$pval, method = 'BH')

#genus vs gpr433
cor.genus.result.gpr43<-list()
for (i in 3:length(cor.genus)){
  genus<-colnames(cor.genus)[i]
  mod <- cor.test(cor.genus[,i], cor.genus$Gpr43, method=c("spearman"))
  pval <- mod$p.value
  r<-mod$estimate
  cor.genus.result.gpr43 <- rbind(cor.genus.result.gpr43, data.frame(genus=genus, 
                                                                     pval=pval,
                                                                     r=r))
}
#adjust for p value by BH method
cor.genus.result.gpr43$FDR_BH <- p.adjust(cor.genus.result.gpr43$pval, method = 'BH')

#cor analysis genus vs scfa
#inport data
cor.genus<-read.csv("taxa_vs_scfa.csv",header=T,row.names = 1)

#spearman cor analysis
cor.genus.scfa.result <-list()
for (j in 1:4) {
    for (i in 5:length(cor.scfa)){
    scfa<-colnames(cor.scfa)[j]
    genus<-colnames(cor.scfa)[i]
    mod <- cor.test(cor.scfa[,i], cor.scfa[,j], method=c("spearman"))
    pval <- mod$p.value
    r<-mod$estimate
   cor.genus.scfa.result <- rbind(cor.genus.scfa.result, data.frame(scfa=scfa,
                                                                    genus=genus, 
                                                                     pval=pval,
                                                                     r=r))
    }
}

#adjust for p value by BH method
cor.genus.scfa.result<- na.omit(cor.genus.scfa.result)  %>% group_by(scfa) %>%  mutate(FDR_BH=p.adjust(pval, method = 'BH')) %>% ungroup()

#rename content and combine tables
C2<-filter(cor.genus.scfa.result,scfa=="C2.Acetate") %>% rename(C2.scfa=scfa,C2.pval=pval,C2=r,C2.FDR_BH=FDR_BH)
C3<-filter(cor.genus.scfa.result,scfa=="C3.Propionate") %>% rename(C3.scfa=scfa,C3.pval=pval,C3=r,C3.FDR_BH=FDR_BH)
C4<-filter(cor.genus.scfa.result,scfa=="C4.Butyrate") %>% rename(C4.scfa=scfa,C4.pval=pval,C4=r,C4.FDR_BH=FDR_BH)
Total_SCFA<-filter(cor.genus.scfa.result,scfa=="Total_SCFA") %>% rename(Total.scfa=scfa,T.pval=pval,Total_SCFA=r,T.FDR_BH=FDR_BH)
SCFA<-merge(C2,C3)%>% merge(.,C4)%>%merge (.,Total_SCFA)
SCFA.pval<-select(SCFA,genus,C2.FDR_BH,C3.FDR_BH,C4.FDR_BH,T.FDR_BH) %>% column_to_rownames(var="genus")
SCFA.star<-ifelse(SCFA.pval<0.01,"**", ifelse (SCFA.pval<0.05, "*","")) %>% t()
SCFA.cor<-select(SCFA,genus,C2,C3,C4,Total_SCFA) %>% column_to_rownames(var="genus") %>% t()

#visualization by pheatmap
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 50, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

paletteLength <- 100
myBreaks <- c(seq(-0.9, 0, length.out=ceiling(paletteLength/2)),seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(paletteLength)

tiff('genus_scfa.tiff', units="in", width=16, height=5.8, res=300, compression = 'lzw')
pheatmap(SCFA.cor,display_numbers = SCFA.star,fontsize =12,fontsize_number=20,
         cluster_rows = F, cluster_cols = F,cellwidth=24,cellheight=24,breaks=myBreaks,color<-myColor,main="Taxa correlated to Fecal SCFA",legend=T)
dev.off()

#cor analysis genus vs gpr41/43
cor.scfa.gprs<-read.csv("scfa_vs_gprs.csv",header=T,row.names = 1)
cor.scfa.gprs.result <-list()
for (j in 1:2) {
  for (i in 3:length(cor.scfa.gprs)){
    gprs<-colnames(cor.scfa.gprs)[j]
    scfa<-colnames(cor.scfa.gprs)[i]
    mod <- cor.test(cor.scfa.gprs[,i], cor.scfa.gprs[,j], method=c("spearman"))
    pval <- mod$p.value
    r<-mod$estimate
    cor.scfa.gprs.result <- rbind(cor.scfa.gprs.result, data.frame(gprs=gprs,
                                                                     scfa=scfa, 
                                                                     pval=pval,
                                                                     r=r))
  }
}

cor.scfa.gprs.result<- cor.scfa.gprs.result  %>% group_by(gprs) %>%  mutate(FDR_BH=p.adjust(pval, method = 'BH')) %>% ungroup()


#cor analysis pathway vs gprs
cor.pw.gprs<-read.csv("pw_vs_gprs.csv",header=T,row.names = 1)

cor.pw.gprs.result <-list()
for (j in 1:2) {
  for (i in 3:length(cor.pw.gprs)){
    gprs<-colnames(cor.pw.gprs)[j]
    pw<-colnames(cor.pw.gprs)[i]
    mod <- cor.test(cor.pw.gprs[,i], cor.pw.gprs[,j], method=c("spearman"))
    pval <- mod$p.value
    r<-mod$estimate
    cor.pw.gprs.result <- rbind(cor.pw.gprs.result, data.frame(gprs=gprs,
                                                                   scfa=scfa, 
                                                                   pval=pval,
                                                                   r=r))
  }
}

cor.pw.gprs.result<- cor.pw.gprs.result  %>% group_by(gprs) %>%  mutate(FDR_BH=p.adjust(pval, method = 'BH')) %>% ungroup()

write.csv(cor.pw.gprs.result, "cor.pw.gprs.result.csv")

cor.carbon<-read.csv("Spearman.carbon.cor.csv",header=T,row.names = 1) %>% t()
padj.carbon<-read.csv("Spearman.carbon.p.csv",header=T,row.names=1)
star.carbon<-ifelse(padj.carbon<0.01,"**", ifelse (padj.carbon<0.05, "*","")) %>% t()

tiff('carbon degradation.tiff', units="in", width=3.2, height=1.6, res=300, compression = 'lzw')
pheatmap(cor.carbon,display_numbers = star.carbon,fontsize =6,fontsize_number=10,cluster_rows = F, cluster_cols = F,cellwidth=12,cellheight=12,breaks=myBreaks,color<-myColor,main="Carbohydrate Degradation",legend=F)
dev.off()

cor.fermentation<-read.csv("Spearman.fermentation.cor.csv",header=T,row.names = 1) %>% t()
padj.fermentation<-read.csv("Spearman.fermentation.p.csv",header=T,row.names=1)
star.fermentation<-ifelse(padj.fermentation<0.01,"**", ifelse (padj.fermentation<0.05, "*","")) %>% t()

tiff('Fermentation.legend.tiff', units="in", width=3.2, height=1.7, res=300, compression = 'lzw')
pheatmap(cor.fermentation, display_numbers = star.fermentation,fontsize =6,fontsize_number=10,cluster_rows = F, cluster_cols = F,cellwidth=12,cellheight=12,breaks=myBreaks,color<-myColor,main="Fermentation",legend_breaks = c(-0.8,-0.4,0,0.4,0.8))
dev.off()
