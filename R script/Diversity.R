#environment
library(tidyverse)
library(tidytext)
library(vegan)
library(phyloseq)
library(ape)
library(philr)
library(data.table)
library(iRF)
library(ROCR)

#combine OTU table
otus=list.files("./OTU",full.names=T)%>%
  lapply(., function(x) suppressMessages(read_tsv(x, skip=1))) %>% reduce(full_join, by = "#OTU ID")
otus[is.na(otus)] <- 0 
otus=column_to_rownames(otus,var="#OTU ID")


#reads across all the samples
p1=data.frame(Counts=colSums(otus)) %>%
  rownames_to_column("SampleID") %>%
  ggplot(aes(x=Counts)) +
  geom_freqpoly() +
  scale_x_log10() +
  theme_bw()+
  theme(panel.background = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5), 
        axis.title= element_text(color="black",size=10), 
        axis.text= element_text(color="black",size=10), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))+
  ylab("# Samples") +
  xlab("Reads") +
  ggtitle("Reads distribution across all samples")

tiff('Fig_final/read across samples.tiff', units="in", width=6, height=5, res=600, compression = 'lzw')
p1
dev.off()


#reads across all the studies
metadata=read.csv("metadata.csv",header=T)
counts<-data.frame(Counts=colSums(otus))%>%rownames_to_column("id")%>%left_join(metadata)
level<-counts %>% group_by(study) %>% summarize(median=median(Counts)) %>% arrange(desc(median)) %>% pull(study)
p2=counts %>% mutate(study=factor(study, levels=level)) %>%
  ggplot(aes(x=study, y=Counts, fill=study)) +
  geom_boxplot(outlier.size=0.1) +
  theme_bw()+ 
  theme(panel.background = element_blank(),
        plot.title = element_text(size=10, hjust = 0.5), 
        axis.title= element_text(color="black",size=10), 
        axis.text= element_text(color="black",size=10), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))+
  scale_y_log10() +
  ylab("Reads") +
  xlab("Study") +
  ggtitle("Read distribution across all studies") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=0, hjust=1))+
  geom_hline(yintercept = 10000, linetype="dashed", color="red")+ theme(axis.text.x=element_text(angle=45, hjust=1))

tiff('Fig_final/read across studies.tiff', units="in", width=6, height=5, res=600, compression = 'lzw')
p2
dev.off()
rm(level, counts)

#import taxonomy
taxonomy<-read_tsv("97_otu_taxonomy.txt", col_names = c("OTU","Taxonomy")) %>%
  separate(Taxonomy, 
           c("Kingdom",
             "Phylum",
             "Class",
             "Order",
             "Family",
             "Genus",
             "Species"
           ), 
           sep="; ", 
           remove=FALSE) %>%
  as.data.frame()
rownames(taxonomy)<-taxonomy$OTU

#generate a list of otu tables with filtering noisy features and subsampled
OTUs<-list()
OTUs$Raw=otus
taxonomy<-subset(taxonomy, OTU %in% rownames(OTUs$Raw))
filtered<-otus[apply(otus,1,function(x){length(grep("TRUE",x!=0))}>=2),] #present min 2 samples
filtered<-filtered[((rowSums(filtered))>=5),] #min 5 reads across the samples
OTUs$Filtered<-filtered
OTUs$Subsampled<-vegan::rrarefy(t(OTUs$Filtered),10000)%>%t() #rarefy to 10K per sample
rm(filtered)

#generate otu tress
trees<-list()
trees$Raw<-read_tree_greengenes("97_otus.tree")
trees$Raw$edge.length[is.nan(trees$Raw$edge.length)]<-0 #convert nan to 0
trees$Raw<-drop.tip(trees$Raw, tip=trees$Raw$tip.label[!trees$Raw$tip.label %in% rownames(OTUs$Raw)])
trees$Subsampled<-drop.tip(trees$Raw, tip=trees$Raw$tip.label[!trees$Raw$tip.label %in% rownames(OTUs$Subsampled)])

#Summarize.Taxa function from jbisanz/MicrobeR
Summarize.Taxa<-function(FEATURES, TAXONOMY){
  
  
  if(sum(!rownames(FEATURES) %in% rownames(TAXONOMY))>0){stop("Mismatch between feature table and taxonomy table detected. There are features which are not in the taxonomy table. Check dimensions and row names.")}#check that all features in feature table are also in taxonomy table
  
  TAXONOMY<-TAXONOMY[rownames(FEATURES),] #order the taxonomy table in the order as the taxonomy table
  
  taxastrings<-data.table(
    FeatureID=rownames(TAXONOMY),
    Kingdom=paste0(TAXONOMY$Kingdom),
    Phylum=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, sep=';'),
    Class=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, sep=';'),
    Order=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, sep=';'),
    Family=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, sep=';'),
    Genus=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, sep=';'),
    Species=paste(TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, TAXONOMY$Species, sep=';'),
    Feature=paste(rownames(TAXONOMY), TAXONOMY$Kingdom, TAXONOMY$Phylum, TAXONOMY$Class, TAXONOMY$Order, TAXONOMY$Family, TAXONOMY$Genus, TAXONOMY$Species, sep=';'),
    stringsAsFactors = F)
  taxastrings[is.na(taxastrings)]<-"NotAssigned"
  
  SUMMARIZED.TAXA=list()
  for (level in colnames(taxastrings)[2:8]){
    SUMMARIZED.TAXA[[level]]<-data.table(FEATURES)
    SUMMARIZED.TAXA[[level]]$TaxaString<-taxastrings[,get(level)] # was important that this was ordered above
    SUMMARIZED.TAXA[[level]]<-as.data.frame(SUMMARIZED.TAXA[[level]][, lapply(.SD, sum), by=TaxaString]) #from data.table vignette, SD is a data.table short hand for the subset of data
    rownames(SUMMARIZED.TAXA[[level]])<-SUMMARIZED.TAXA[[level]]$TaxaString
    SUMMARIZED.TAXA[[level]]$TaxaString<-NULL
  }
  SUMMARIZED.TAXA$Feature<-FEATURES
  rownames(SUMMARIZED.TAXA$Feature)<-taxastrings$Feature
  
  if(diff(range(lapply(SUMMARIZED.TAXA, sum)))!=0){stop("Reads went missing during summary. Debugging will be necesssary.")}
  
  return(SUMMARIZED.TAXA)
}
OTUs$TaxaSummary<-Summarize.Taxa(OTUs$Filtered, taxonomy)

#calculate relative abundance at phylum level
phl<-OTUs$TaxaSummary$Phylum
for(col in names(phl)) {
  phl[paste0(col)] = (phl[col] / sum(phl[col])*100)
  rm(col)}

phab<-
  as.data.frame(phl)%>% rownames_to_column("Phylum") %>%  
  mutate(Phylum=gsub("..+p__","", Phylum)) %>% 
  gather(-Phylum, key="id", value="Abundance")

tp<-
  phab %>% group_by(Phylum) %>% summarize(mean=mean(Abundance)) %>%
  arrange(desc(mean)) %>%
  top_n(7, mean) %>%
  bind_rows(., tibble(Phylum="Other", mean=0))

por<-
  phab %>%
  filter(Phylum=="Bacteroidetes") %>%
  arrange(desc(Abundance)) %>%
  pull(id)

p3=phab %>%
  mutate(Phylum=if_else(Phylum %in% tp$Phylum, Phylum, "Other")) %>%
  mutate(Phylum=factor(Phylum, levels = rev(tp$Phylum))) %>%
  left_join(metadata) %>%
  mutate(Sample=factor(id, levels=por)) %>% mutate(age=factor(age,level=c("young","old"))) %>%
  ggplot(aes(x=Sample, y=Abundance, fill=Phylum, width=1)) +
  geom_bar(stat="identity") +
  facet_grid(~sex+as.factor(age), scales="free", space = "free")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank()) +
  theme(axis.text.y = element_text(color="black",size=8),axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y=element_text(size=8), 
        strip.text = element_text(size=8))+
  scale_fill_manual(values=rev(c(
    "blue4",
    "olivedrab",
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "grey")))

tiff('Fig_final/phylumbar1.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')
p3
dev.off()

rm(phab, phl, tp, por)

#Distance Matrix Generation
DistanceMatrices<-list()
DistanceMatrices$Subsampled<-list()

prop=OTUs$Subsampled
prop=as.data.frame(prop)
for(col in names(prop)) {
  prop[paste0(col)] = (prop[col] / sum(prop[col])*100)
  rm(col)}


DistanceMatrices$Subsampled[["weighted UniFrac"]]<-
  UniFrac(phyloseq(otu_table(prop, taxa_are_rows = T), phy=trees$Subsampled), weighted=T, parallel=T)
DistanceMatrices$Subsampled[["unweighted UniFrac"]]<-
  UniFrac(phyloseq(otu_table(prop, taxa_are_rows = T), phy=trees$Subsampled), weighted=F, parallel=T)

rm(prop)

#PCoA
mousesamps<-metadata %>%  pull(id)
PCoAs<-list()
for (y in names(DistanceMatrices)){
  for(x in names(DistanceMatrices[[y]])){
    tm<-mousesamps[mousesamps %in% colnames(as.matrix(DistanceMatrices[[y]][[x]]))]
    message(x,"-", y)
    PCoAs[[paste0(x,"-",y)]]<-
      as.matrix(DistanceMatrices[[y]][[x]])[tm,tm] %>%
      as.dist(.) %>%
      ape::pcoa()
  }
}
rm(tm, y, x)

plotorder<-c("unweighted UniFrac","weighted UniFrac")

pc=lapply(names(PCoAs), function(x) PCoAs[[x]]$vectors %>% 
            as.data.frame() %>% 
            rownames_to_column("id") %>% 
            dplyr::select(id, Axis.1, Axis.2, Axis.3) %>% 
            mutate(Metric=x)) %>%
  do.call(bind_rows, .) %>%
  as_tibble() %>%
  left_join(metadata) %>%
  separate(Metric, c("Metric","Dataset"), sep="-") %>%
  mutate(Metric=factor(Metric, levels=plotorder)) 
rm(plotorder)

PCvar<-lapply(names(PCoAs), function(x){
  PCoAs[[x]]$values %>%
    as.data.frame() %>%
    rownames_to_column("Axis") %>%
    as_tibble() %>%
    mutate(Metric=x) %>%
    mutate(Axis=as.numeric(Axis)) %>%
    filter(Axis<=10)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(Pvar=round (100*Relative_eig,2)) %>%
  separate(Metric, c("Metric","Dataset"), sep="-")


p4=ggplot(pc,aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=age,shape=study), size=3) +
  facet_wrap(~Metric, scales="free", ncol=2) +theme_minimal()+
  theme(panel.border = element_blank(), axis.line = element_line(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("") +
  ylab("") + scale_color_manual(values=c("firebrick3","dodgerblue3"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2), 
        panel.background = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5), 
        axis.title= element_text(color="black",size=10), 
        axis.text= element_text(color="black",size=10), 
        panel.border = element_rect(colour = "black", fill=NA))

PC1_unw <- PCvar %>% filter(., Metric=="unweighted UniFrac", Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- PCvar %>% filter(., Metric=="unweighted UniFrac", Axis=="2") %>% select(Pvar) %>% as.character()

p4_unweight<-ggplot(filter(pc,Metric=="unweighted UniFrac"),aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=age,shape=study), size=3) +theme_minimal()+
  theme(panel.border = element_blank(), axis.line = element_line(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) +
  ylab(paste("PCoA2"," (", paste(PC2_unw), "%)")) + scale_color_manual(values=c("firebrick3","dodgerblue3"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2), 
        panel.background = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5), 
        axis.title= element_text(color="black",size=10), 
        axis.text= element_text(color="black",size=10), 
        panel.border = element_rect(colour = "black", fill=NA))

tiff('Fig_final/PCoA unweighted.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')
p4_unweight+ggtitle("unweighted UniFrac")
dev.off()

PC1_w <- PCvar %>% filter(., Metric=="weighted UniFrac", Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- PCvar %>% filter(., Metric=="weighted UniFrac", Axis=="2") %>% select(Pvar) %>% as.character()

p4_weight<-ggplot(filter(pc,Metric=="weighted UniFrac"),aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(color=age,shape=study), size=3) +theme_minimal()+
  theme(panel.border = element_blank(), axis.line = element_line(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) +
  ylab(paste("PCoA2"," (", paste(PC2_w), "%)")) + scale_color_manual(values=c("firebrick3","dodgerblue3"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2), 
        panel.background = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5), 
        axis.title= element_text(color="black",size=10), 
        axis.text= element_text(color="black",size=10), 
        panel.border = element_rect(colour = "black", fill=NA))

tiff('Fig_final/PCoA weighted.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')
p4_weight+ggtitle("weighted UniFrac")
dev.off()

rm(PC1_unw,PC2_unw,PC1_w,PC2_w,PCvar,pc)

#ADONIS
ADONIS<-tibble()
for(x in names(DistanceMatrices)){
  for(y in names(DistanceMatrices[[x]])){
    ts<-mousesamps[mousesamps %in% labels(DistanceMatrices[[x]][[y]])]
    tm<-metadata %>% as.data.frame() %>% remove_rownames() %>% filter(id %in% ts) %>% column_to_rownames("id") 
    td<-as.matrix(DistanceMatrices[[x]][[y]])[ts,ts] %>% as.dist()
    tm<-tm[labels(td),]
    ta<-adonis(td~age, data=tm, strata=tm$study, permutations=999, parallel=1)
    ta<-ta$aov.tab %>% 
      as.data.frame() %>% 
      rownames_to_column("Term") %>%
      mutate(Metric=paste0(x,"-",y)) %>% 
      dplyr::select(Metric, everything()) %>%
      rename(Pvalue=`Pr(>F)`)
    ADONIS<-bind_rows(ADONIS, ta)
    rm(td, tm, ta, ts)
    gc()
  }
}

rm(x,y)

#perstudy otu table
processStudy<-function(Study){
  PS<-list()
  PS$Study<-Study
  PS$Metadata<-subset(metadata, study==Study)
  PS$OTUs_raw<-OTUs$Raw[,PS$Metadata$id]
  PS$OTUs_raw<-PS$OTUs_raw[rowSums(PS$OTUs_raw)>=1,]
}

PerStudy<-list()
for(study in (unique(metadata$study))){
  PerStudy[[study]]<-processStudy(study)
}

#perstudy phylum FB ratio
per.taxa<-list()
for(study in (unique(metadata$study))){
  per.taxa[[study]]<-Summarize.Taxa(PerStudy[[study]],taxonomy %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("OTU"))
}

per.phylum<-list()
for(study in (unique(metadata$study))){
  per.phylum[[study]]<-per.taxa[[study]]$Phylum
  for(col in names(per.phylum[[study]])) {
    per.phylum[[study]][paste0(col)] = (per.phylum[[study]][col] / sum(per.phylum[[study]][col])*100)
    rm(col)}
  per.phylum[[study]]<-as.data.frame(t(per.phylum[[study]]))
  per.phylum[[study]]$'FB ratio'=per.phylum[[study]]$`k__Bacteria;p__Firmicutes`/per.phylum[[study]]$`k__Bacteria;p__Bacteroidetes`
  rm(study)
}

per.fb<-list()
for(study in (unique(metadata$study))){
  per.fb[[study]]<-per.phylum[[study]] %>% as.data.frame() %>% rownames_to_column("id") %>% dplyr::select(id, 'FB ratio')
}


per.fb.utest<-list()
for(study in (unique(metadata$study))){
  per.fb.utest[[study]] <- per.fb[[study]] %>% left_join(metadata) %>%
    do(broom::tidy(wilcox.test(`FB ratio`~age, data=., conf.int=TRUE, conf.level=0.95))
    ) %>%
    mutate(study=study) %>%
    dplyr::select(study, `mean of difference`=estimate, Pvalue=p.value, 
                  CI_low=conf.low, CI_high=conf.high)
}


fb<-do.call(bind_rows,per.fb) %>% left_join(metadata) %>% mutate(age=factor(age, levels=c("young","old"))) 
fit.fb<-lmerTest:::lmer(`FB ratio`~age+(1|study), data=fb)
cf.fb.all<-confint(fit.fb)
fb.comb.all<-tibble(study=character(0), `mean of difference`=numeric(0), 
                    Pvalue=numeric(0), CI_low=numeric(0), CI_high=numeric(0))
fb.comb.all<-tibble(study="Combined (All, n=134)", 
                    `mean of difference`=summary(fit.fb)$coefficients["ageold", "Estimate"],
                    Pvalue=anova(fit.fb)$`Pr(>F)`,CI_low=cf.fb.all["ageold",1], 
                    CI_high=cf.fb.all["ageold",2])

p5<-do.call(bind_rows,per.fb.utest) %>% bind_rows(fb.comb.all)%>%
  mutate(Significant=case_when(
    Pvalue<0.05 & `mean of difference`>0 ~ "increased in old",
    Pvalue<0.05 & `mean of difference`<0 ~ "increased in young",
    TRUE~"not significant"
  ))%>%
  ungroup() %>%
  ggplot(aes(x=`mean of difference`, y=study, color=Significant)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
  geom_point() +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1),axis.text = element_text(color="black"))+
  scale_color_manual(values=c("firebrick3","dodgerblue3","black"))+xlab("mean of difference (95% CI)")

tiff('Fig_final/FB Forest.tiff', units="in", width=4.5, height=4.5, res=600, compression = 'lzw')
p5+ggtitle("FB ratio")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

rm(fb,per.taxa,per.fb,per.phylum,fb.allstudy.comb,fit.fb,cf.fb.all,fb.allstudy.comb)

#perstudy filtered otu table and tree
#filtering: present min 1 samples & min 3 reads across the samples per study
per.filter<-list()
for(study in (unique(metadata$study))){
  per.filter[[study]]<-PerStudy[[study]][apply(PerStudy[[study]],1,function(x){length(grep("TRUE",x!=0))}>=1),]
  per.filter[[study]]<-per.filter[[study]][((rowSums(per.filter[[study]]))>=3),]
  rm(study)
}

per.depth<-list()
for(study in (unique(metadata$study))){
  per.depth[[study]]=min(colSums(per.filter[[study]]))
  rm(study)
}

per.subsample<-list()
for(study in (unique(metadata$study))){
  per.subsample[[study]]<-vegan::rrarefy(t(per.filter[[study]]),per.depth[[study]])%>%t()
}

per.tree<-list()
for(study in (unique(metadata$study))){
  per.tree[[study]]<-drop.tip(trees$Raw, trees$Raw$tip.label[!trees$Raw$tip.label %in% rownames(per.subsample[[study]])])
}

#perstudy alpha diversity
per.AD<-list()
for(study in (unique(metadata$study))){
  per.AD[[study]]<-data.frame(
    Shannon=vegan::diversity(per.subsample[[study]],index="shannon",MARGIN = 2),
    Simpson=vegan::diversity(per.subsample[[study]],index="simpson",MARGIN = 2),
    Chao1=vegan::estimateR(t(per.subsample[[study]])) %>% t() %>% as.data.frame %>% rename(Chao1=S.chao1) %>% pull(Chao1),
    FaithsPD=picante::pd(t(per.subsample[[study]]), per.tree[[study]], include.root=F)$PD)
}


per.AD.stats<-list()
for(study in (unique(metadata$study))){
  per.AD.stats[[study]]<-per.AD[[study]] %>% rownames_to_column("id") %>%
    dplyr::select(id, everything())%>%
    as_tibble() %>%
    gather(-id, key="Metric", value="Diversity") %>%
    left_join(metadata[,c("id","age", "study","sex")]) %>%
    dplyr::select(id, study, age, sex,Metric, Diversity)
}

per.AD.stats.log<-list()
for(study in (unique(metadata$study))){
  per.AD.stats.log[[study]]<-per.AD.stats[[study]] %>%
    left_join(per.AD.stats[[study]]%>%
                group_by(Metric, age) %>%
                summarize(mean=mean(log2(Diversity))) %>%
                spread(key=age, value=mean) %>%
                rename(mean_log2old=old, mean_log2young=young)) %>%
    mutate(log2FC=log2(Diversity)-mean_log2young)
}


per.AD.stats.ttest<-list()
for(study in (unique(metadata$study))){
  per.AD.stats.ttest[[study]] <- per.AD.stats.log[[study]] %>%
    group_by(Metric) %>% 
    do(
      broom::tidy(t.test(log2FC~age, data=., conf.int=TRUE, conf.level=0.95))
    ) %>%
    mutate(study=study) %>%
    dplyr::select(study, Metric, log2FC=estimate, Pvalue=p.value, 
                  mean_young=estimate2, mean_old=estimate1, 
                  CI_low=conf.low, CI_high=conf.high)
}

alpha.all<-do.call(bind_rows,per.AD.stats.log) %>% mutate(age=factor(age, levels=c("young","old")))

comb.all<-tibble(study=character(0), Metric=character(0), log2FC=numeric(0), Pvalue=numeric(0), mean_young=numeric(0), mean_old=numeric(0), CI_low=numeric(0), CI_high=numeric(0))
for(i in unique(alpha.all$Metric)){
  fit<-lmerTest:::lmer(log2FC~age+(1|study), data=subset(alpha.all, Metric==i))
  cf<-confint(fit,level = 0.95)
  comb.all<-bind_rows(comb.all, tibble(
    study="Combined (All, n=134)", 
    Metric=i, 
    log2FC=summary(fit)$coefficients["ageold", "Estimate"], 
    Pvalue=anova(fit)$`Pr(>F)`, 
    mean_young=NA, 
    mean_old=NA, 
    CI_low=cf["ageold",1], 
    CI_high=cf["ageold",2]
  ))
}


all.AD<-do.call(bind_rows,per.AD.stats.ttest)
p6<-bind_rows(all.AD,comb.all)%>%
  mutate(Significant=case_when(
    Pvalue<0.05 & log2FC>0 ~ "increased in old",
    Pvalue<0.05 & log2FC<0 ~ "young",
    TRUE~"not significant"
  ))%>%
  ungroup() %>%
  mutate(Metric=factor(Metric, 
                       levels=c("Shannon","Simpson", "FaithsPD","Chao1"))) %>%
  ggplot(aes(x=log2FC, y=study, color=Significant)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
  geom_point() +
  facet_grid(~Metric, scales="free_x") +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_color_manual(values=c("firebrick3","black","dodgerblue3"))+xlab("log2FC (95% CI)")

tiff('Fig_final/Alpha Forest.tiff', units="in", width=8, height=4.5, res=600, compression = 'lzw')
p6
dev.off()

#perstudy beta diversity
per.otu.prop<-per.subsample
for(study in (unique(metadata$study))){
  per.otu.prop[[study]]<-as.data.frame(per.otu.prop[[study]])
}

for(study in (unique(metadata$study))){
  for(col in names(per.otu.prop[[study]])) {
    per.otu.prop[[study]][paste0(col)] = (per.otu.prop[[study]][col] / sum(per.otu.prop[[study]][col])*100)
    rm(col)}
  per.otu.prop[[study]]=t(per.otu.prop[[study]])
}


#perstudy PCoA analysis
per.BD<-list()  
for(study in (unique(metadata$study))){
  per.BD[[study]]<-list()
  per.BD[[study]][["unweighted UniFrac"]]<-
    UniFrac(phyloseq(otu_table(t(per.otu.prop[[study]]), taxa_are_rows = T), phy=per.tree[[study]]), weighted=F, parallel=T)
  per.BD[[study]][["weighted UniFrac"]]<-
    UniFrac(phyloseq(otu_table(t(per.otu.prop[[study]]), taxa_are_rows = T), phy=per.tree[[study]]), weighted=T, parallel=T)
}  

mousesamps<-metadata %>%  pull(id)

per.pcoa<-list()
for(study in (unique(metadata$study))){
  for (x in names(per.BD[[study]])){
    tm<-mousesamps[mousesamps %in% colnames(as.matrix(per.BD[[study]][[x]]))]
    per.pcoa[[study]][[paste0(x)]]<-
      as.matrix(per.BD[[study]][[x]])[tm,tm] %>%
      as.dist(.) %>%
      ape::pcoa()
  }
}


#unweighted PCoA plot
per.pcoa.plot_unw<-list()
for(study in (unique(metadata$study))){
  per.pcoa.plot_unw[[study]]<-lapply(names(per.pcoa[[study]]), function(x) per.pcoa[[study]][[x]]$vectors %>% 
                                       as.data.frame() %>% 
                                       rownames_to_column("id") %>% 
                                       dplyr::select(id, Axis.1, Axis.2, Axis.3) %>% 
                                       mutate(Metric=x)) %>%
    do.call(bind_rows, .) %>%
    as_tibble() %>%
    left_join(metadata) %>% filter(Metric=="unweighted UniFrac") %>%
    ggplot(aes(x=Axis.1, y=Axis.2)) +
    geom_point(aes(color=age), shape=16, alpha=0.8) +theme_minimal()+
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") + scale_color_manual(values=c("firebrick3","dodgerblue3"))+
    theme(legend.key = element_rect(fill = NA, colour = NA, size = 0.25), 
          panel.background = element_blank(),
          plot.title = element_text(size=10, hjust = 0.5), 
          axis.title= element_text(color="black",size=10), 
          axis.text= element_text(color="black",size=10), 
          panel.border = element_rect(colour = "black", fill=NA, size=1.5))
}


per.pcoa.pc_unw<-list()
for(study in (unique(metadata$study))){
  per.pcoa.pc_unw[[study]][["unweighted UniFrac"]]<-per.pcoa[[study]][["unweighted UniFrac"]]$values%>%
    as.data.frame() %>%
    rownames_to_column("Axis") %>%
    as_tibble() %>%
    mutate(Metric="unweighted UniFrac") %>%
    mutate(Axis=as.numeric(Axis)) %>%
    filter(Axis<=8)%>%
    mutate(Pvar=round(100*Relative_eig,2)) %>%
    mutate(study=study)
}


#Exp1 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp1 (Female, n=10)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp1 (Female, n=10)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp1.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp1 (Female, n=10)`+ ggtitle("Exp1 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()


#Exp2 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp2 (Female, n=12)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp2 (Female, n=12)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp2.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp2 (Female, n=12)`+ ggtitle("Exp2 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()

#Exp3 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp3 (Female, n=10)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp3 (Female, n=10)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp3.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp3 (Female, n=10)`+ ggtitle("Exp3 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()

#Exp4 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp4 (Male, n=16)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp4 (Male, n=16)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp4.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp4 (Male, n=16)`+ ggtitle("Exp4 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()

#Exp5 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp5 (Male, n=65)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp5 (Male, n=65)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp5.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp5 (Male, n=65)`+ ggtitle("Exp5 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()

#Exp6 PCoA Unweighted
PC1_unw <- bind_rows(per.pcoa.pc_unw$`Exp6 (Male, n=21)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_unw <- bind_rows(per.pcoa.pc_unw$`Exp6 (Male, n=21)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp6.unw.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_unw$`Exp6 (Male, n=21)`+ ggtitle("Exp6 unweighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_unw), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_unw), "%)"))
dev.off()


#weighted PCoA plot
per.pcoa.plot_w<-list()
for(study in (unique(metadata$study))){
  per.pcoa.plot_w[[study]]<-lapply(names(per.pcoa[[study]]), function(x) per.pcoa[[study]][[x]]$vectors %>% 
                                     as.data.frame() %>% 
                                     rownames_to_column("id") %>% 
                                     dplyr::select(id, Axis.1, Axis.2, Axis.3) %>% 
                                     mutate(Metric=x)) %>%
    do.call(bind_rows, .) %>%
    as_tibble() %>%
    left_join(metadata) %>% filter(Metric=="weighted UniFrac") %>%
    ggplot(aes(x=Axis.1, y=Axis.2)) +
    geom_point(aes(color=age), shape=16, alpha=0.8) +theme_minimal()+
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") + scale_color_manual(values=c("firebrick3","dodgerblue3"))+
    theme(legend.key = element_rect(fill = NA, colour = NA, size = 0.25), 
          panel.background = element_blank(),
          plot.title = element_text(size=10, hjust = 0.5), 
          axis.title= element_text(color="black",size=10), 
          axis.text= element_text(color="black",size=10), 
          panel.border = element_rect(colour = "black", fill=NA, size=1.5))
}


per.pcoa.pc_w<-list()
for(study in (unique(metadata$study))){
  per.pcoa.pc_w[[study]][["weighted UniFrac"]]<-per.pcoa[[study]][["weighted UniFrac"]]$values%>%
    as.data.frame() %>%
    rownames_to_column("Axis") %>%
    as_tibble() %>%
    mutate(Metric="weighted UniFrac") %>%
    mutate(Axis=as.numeric(Axis)) %>%
    filter(Axis<=8)%>%
    mutate(Pvar=round(100*Relative_eig,2)) %>%
    mutate(study=study)
}

#Exp1 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp1 (Female, n=10)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp1 (Female, n=10)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp1.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp1 (Female, n=10)`+ ggtitle("Exp1 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()


#Exp2 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp2 (Female, n=12)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp2 (Female, n=12)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp2.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp2 (Female, n=12)`+ ggtitle("Exp2 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()

#Exp3 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp3 (Female, n=10)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp3 (Female, n=10)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp3.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp3 (Female, n=10)`+ ggtitle("Exp3 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()

#Exp4 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp4 (Male, n=16)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp4 (Male, n=16)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp4.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp4 (Male, n=16)`+ ggtitle("Exp4 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()

#Exp5 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp5 (Male, n=65)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp5 (Male, n=65)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp5.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp5 (Male, n=65)`+ ggtitle("Exp5 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()

#Exp6 PCoA weighted
PC1_w <- bind_rows(per.pcoa.pc_w$`Exp6 (Male, n=21)`) %>% filter(., Axis=="1") %>% select(Pvar) %>% as.character()
PC2_w <- bind_rows(per.pcoa.pc_w$`Exp6 (Male, n=21)`) %>% filter(., Axis=="2") %>% select(Pvar) %>% as.character()
tiff('Fig_final/Exp6.weighted.tiff', units="in", width=4, height=3, res=600, compression = 'lzw')
per.pcoa.plot_w$`Exp6 (Male, n=21)`+ ggtitle("Exp6 weighted UniFrac")+ xlab(paste("PCoA1"," (", paste(PC1_w), "%)")) + ylab(paste("PCoA2"," (", paste(PC2_w), "%)"))
dev.off()

#perstudy adonis
per.ADONIS<-list()
for(study in (unique(metadata$study))){
  per.ADONIS[[study]]<-tibble()
  for(y in names(per.BD[[study]])){
    ts<-mousesamps[mousesamps %in% labels(per.BD[[study]][[y]])]
    tm<-metadata %>% as.data.frame() %>% remove_rownames() %>% filter(id %in% ts) %>% column_to_rownames("id") 
    td<-as.matrix(per.BD[[study]][[y]])[ts,ts] %>% as.dist()
    tm<-tm[labels(td),]
    ta<-adonis(td~age, data=tm, permutations=999, parallel=1)
    ta<-ta$aov.tab %>% 
      as.data.frame() %>% 
      rownames_to_column("Term") %>%
      mutate(Metric=paste0(y)) %>% 
      dplyr::select(Metric, everything()) %>%
      rename(Pvalue=`Pr(>F)`)
    per.ADONIS[[study]]<-bind_rows(per.ADONIS[[study]], ta)
    gc()
  }
}
rm(td, tm, ta, ts)

adonis.per.study<-bind_rows(per.ADONIS, .id="Study")








