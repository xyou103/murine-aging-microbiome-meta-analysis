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

#Data partition training and testing data
set.seed(123)
metadata<-read.csv("metadata.csv",header=T)
RF.metadata<-metadata %>% 
  mutate(age=factor(age, levels=c("old","young"))) %>% as.data.frame()
rownames(RF.metadata)<-RF.metadata$id

TrainingSet<-RF.metadata %>% pull(id) %>% sample (.,round(2/3*length(.)),0)

RF.metadata<-
  RF.metadata %>%
  mutate(Task=case_when(
    id %in% TrainingSet ~ "Training",
    TRUE ~ "Testing"
  ))

RF.b6.metadata<-read.csv("B6-metadata.csv",header=T)
RF.exp7.metadata<-read.csv("Exp7-metadata.csv",header=T)
RF.metadata=merge(RF.metadata,RF.b6.metadata,all=T)
RF.metadata<-merge(RF.metadata,RF.exp7.metadata,all=T)

RF.metadata %>%
  group_by(Task) %>%
  summarize(Nsample=length(id)) 

#OTU table import for external data B6 metagenome
B6_otu<-read_tsv("B6-meta-table.txt",skip=1) %>% column_to_rownames("#OTU ID")
B6_filtered<-B6_otu[apply(B6_otu,1,function(x){length(grep("TRUE",x!=0))}>=2),] #present min 2 samples
B6_filtered<-B6_filtered[((rowSums(B6_filtered))>=5),] #min 5 reads across the samples
B6_taxa<-Summarize.Taxa(B6_filtered,taxonomy)
B6_RA<-B6_filtered
for(col in names(B6_RA)) {
    B6_RA[paste0(col)] = (B6_RA[col]/sum(B6_RA[col]))
    rm(col)}
B6_RA<-as.data.frame(B6_RA) %>% rownames_to_column ("taxa")

#OTU table import external data for Exp7 
Exp7_otu<-read_tsv("Exp7-table.txt",skip=1) %>% column_to_rownames("#OTU ID")
Exp7_filtered<-Exp7_otu[apply(Exp7_otu,1,function(x){length(grep("TRUE",x!=0))}>=2),] #present min 2 samples
Exp7_filtered<-Exp7_filtered[((rowSums(Exp7_filtered))>=5),] #min 5 reads across the samples
Exp7_taxa<-Summarize.Taxa(Exp7_filtered,taxonomy)
Exp7_RA<-Exp7_filtered
for(col in names(Exp7_RA)) {
  Exp7_RA[paste0(col)] = (Exp7_RA[col]/sum(Exp7_RA[col]))
  rm(col)}
Exp7_RA<-as.data.frame(Exp7_RA) %>% rownames_to_column ("taxa")

#OTU random forest
#OTUs a list is defined in file Diversity.R 
OTU.RA<-OTUs$Filtered %>% as.data.frame()
for(col in names(OTU.RA)) {
  OTU.RA[paste0(col)] = (OTU.RA[col]/sum(OTU.RA[col]))
  rm(col)}
OTU.RA<-OTU.RA %>% rownames_to_column ("taxa")
OTU.RA.full<-inner_join(OTU.RA,B6_RA,by="taxa") %>% inner_join(.,Exp7_RA, by="taxa") %>% 
  column_to_rownames("taxa")

RF.otu.RA<-iRF::iRF(x=t(OTU.RA.full[,subset(RF.metadata, Task=="Training")$id]),
                     y=subset(RF.metadata, Task=="Training")$age, ntree=500,n.iter =10)

RF.otu.RA$rf.list #check OOB error rate and select the smallest one 


Predictions<-list()
AUCs<-list()

for(Tasks in unique(RF.metadata$Task)){
  tm<-predict(RF.otu.RA$rf.list[[10]], newdata=t(OTU.RA.full[,subset(RF.metadata, Task==Tasks)$id]), type="prob")[,2] %>% 
    prediction(., subset(RF.metadata, Task==Tasks)$age) %>%
    performance(., "tpr","fpr")
  Predictions[[Tasks]]<-
    tibble(Task=Tasks, FPR=unlist(tm@x.values), TPR=unlist(tm@y.values))
  
  auc<-predict(RF.otu.RA$rf.list[[10]], newdata=t(OTU.RA.full[,subset(RF.metadata, Task==Tasks)$id]), type="prob")[,2] %>%
    prediction(., subset(RF.metadata, Task==Tasks)$age) %>%
    performance(., "auc")
  
  AUCs[[Tasks]]<-
    tibble(
      Task=Tasks,
      AUC=auc@y.values[[1]]
    )
}

Predictions<-lapply(Predictions, function(x) do.call(bind_rows, x)) %>% do.call(bind_rows, .)
AUCs<-lapply(AUCs, function(x) do.call(bind_rows, x)) %>% do.call(bind_rows, .)

roc.ra<-Predictions %>% mutate(Task = replace(Task, Task == "Testing", "Test")) %>%
  mutate(Dataset=factor(Task, levels=c("Training","Test","CB6F1", "B6-metagenomics"))) %>%
  ggplot(aes(x=FPR, y=TPR, color=Dataset, group=Dataset)) +
  geom_abline(color="grey", linetype="dashed") + theme_minimal()+ theme(axis.line.x = element_line(size = 1.2, color="grey"))+
  geom_line(size=1.2) +
  xlab("False positive") +
  ylab("True positive") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  coord_cartesian(expand=F)+ggtitle("ROC-OTU")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18, color="black"), axis.text=element_text(size=18,color="black"), legend.text = element_text(size=18))

auc.ra<-AUCs %>% mutate(Task = replace(Task, Task == "Testing", "Test")) %>%
  mutate(Task=factor(Task, levels=c("Training","Test","CB6F1", "B6-metagenomics"))) %>%
  ggplot(aes(x=Task, y=AUC)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  coord_cartesian(ylim=c(0,1), expand=0) +
  scale_fill_viridis_d()+theme_minimal()+ggtitle("AUC-OTU")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18, color="black"), axis.text=element_text(size=18, color="black"),
        axis.text.x = element_text(angle = 30,hjust=1),axis.title.x = element_blank())

tiff("Fig_final/RF.ROC.ra.tiff", units = "in", height=6, width = 6.5, res=300, compression = "lzw")
roc.ra
dev.off()

tiff("Fig_final/RF.AUC.ra.tiff", units = "in", height=6, width = 4.5, res=300, compression = "lzw")
auc.ra
dev.off()

OTU.logFC<-OTU.RA.full %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  as_tibble() %>%
  gather(-OTU, key=id, value=Abundance) %>%
  left_join(RF.metadata %>% dplyr::select(id, Task, age)) %>%
  group_by(OTU, age) %>%
  summarize(mean=mean(Abundance)) %>%
  spread(key=age, value=mean) %>%
  mutate(logFC=log2(old)-log2(young)) %>%
  ungroup() 

RF.features.otu.full<-RF.otu.RA$rf.list[[10]]$importance %>% as.data.frame() %>% arrange(desc(MeanDecreaseGini)) %>% 
  top_n(20, MeanDecreaseGini) %>% rownames_to_column("OTU") %>% left_join(taxonomy, by='OTU') %>% left_join(OTU.logFC, by='OTU')

tiff("Fig_final/RF.feature.otu.tiff", units = "in", height=8, width = 10, res=300, compression = "lzw")
RF.features.otu.full %>% arrange(desc(MeanDecreaseGini)) %>% 
  top_n(20, MeanDecreaseGini) %>% ggplot(aes(x=OTU,y=MeanDecreaseGini,fill=logFC))+
  geom_bar(stat="identity") +theme_minimal() + coord_flip()+ scale_fill_gradient2(low="dodgerblue3", high="firebrick3",n.breaks=8) + 
  xlab("Pathway")+theme(text = element_text(size = 15), axis.text = element_text(color = "black"),axis.title.y = element_blank())+
  ggtitle("Age-specific OTU")
dev.off()

write.csv(RF.features.otu.full,"OTU.feature.csv")

#PW random forest
PW=list.files("./Fun/PW",full.names=T) %>% lapply(., function(x) suppressMessages(read_tsv(x) %>% dplyr::select(,-description))) %>% reduce(full_join,by = "pathway") 
PW[is.na(PW)] <- 0 
PW=column_to_rownames(PW,var="pathway")
PW.filter<-PW[apply(PW,1,function(x){length(grep("TRUE",x!=0))}>=5),]

PW.RA<-PW.filter
for(col in names(PW.RA)) {
  PW.RA[paste0(col)] = (PW.RA[col]/sum(PW.RA[col]))
  rm(col)}

PW.RA<-PW.RA %>% rownames_to_column("pathway")
PW.b6<-read_tsv("./Fun/B6_path_unstrat_des.tsv") %>% dplyr::select(,-description)
PW.RA.b6<-inner_join(PW.RA,PW.b6,by="pathway")
PW.RA.b6[is.na(PW.RA.b6)] <- 0 
PW.exp7<-read_tsv("./Fun/Exp7-PW.tsv") %>% dplyr::select(,-description)
PW.RA.exp7<-PW.exp7
PW.RA.exp7=column_to_rownames(PW.RA.exp7,var="pathway")
for(col in names(PW.RA.exp7)) {
  PW.RA.exp7[paste0(col)] = (PW.RA.exp7[col]/sum(PW.RA.exp7[col]))
  rm(col)}
PW.RA.exp7<-PW.RA.exp7 %>% rownames_to_column("pathway")
PW.full<-inner_join(PW.RA.b6,PW.RA.exp7,by="pathway")

PW.full<-column_to_rownames(PW.full,var="pathway")

PW.des<-list.files("./Fun/PW",full.names=T) %>%
  lapply(., function(x) suppressMessages(read_tsv(x) %>% dplyr::select(feature="pathway",description))) %>% 
  reduce(full_join,by = c("feature","description")) %>% mutate(Type="pathway")

RF.pw.full<-iRF::iRF(x=t(PW.full[,subset(RF.metadata, Task=="Training")$id]),
                     y=subset(RF.metadata, Task=="Training")$age, ntree=500,n.iter =10)

RF.pw.full$rf.list #check OOB rate and select the smallest

Predictions<-list()
AUCs<-list()

for(Tasks in unique(RF.metadata$Task)){
  tm<-predict(RF.pw.full$rf.list[[10]], newdata=t(PW.full[,subset(RF.metadata, Task==Tasks)$id]), type="prob")[,2] %>% 
    prediction(., subset(RF.metadata, Task==Tasks)$age) %>%
    performance(., "tpr","fpr")
  Predictions[[Tasks]]<-
    tibble(Task=Tasks, FPR=unlist(tm@x.values), TPR=unlist(tm@y.values))
  
  auc<-predict(RF.pw.full$rf.list[[10]], newdata=t(PW.full[,subset(RF.metadata, Task==Tasks)$id]), type="prob")[,2] %>%
    prediction(., subset(RF.metadata, Task==Tasks)$age) %>%
    performance(., "auc")
  
  AUCs[[Tasks]]<-
    tibble(
      Task=Tasks,
      AUC=auc@y.values[[1]]
    )
}

Predictions<-lapply(Predictions, function(x) do.call(bind_rows, x)) %>% do.call(bind_rows, .)
AUCs<-lapply(AUCs, function(x) do.call(bind_rows, x)) %>% do.call(bind_rows, .)

roc.pw.b6<-Predictions %>% mutate(Task = replace(Task, Task == "Testing", "Test")) %>%
  mutate(Dataset=factor(Task, levels=c("Training","Test","CB6F1", "B6-metagenomics"))) %>%
  ggplot(aes(x=FPR, y=TPR, color=Dataset, group=Dataset)) +
  geom_abline(color="grey", linetype="dashed") + theme_minimal()+ theme(axis.line.x = element_line(size = 1.2, color="grey"))+
  geom_line(size=1.2) +
  xlab("False positive") +
  ylab("True positive") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  coord_cartesian(expand=F)+ggtitle("ROC-Pathway")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18, color="black"), axis.text=element_text(size=18,color="black"), legend.text = element_text(size=18))

auc.pw.b6<-AUCs %>% mutate(Task = replace(Task, Task == "Testing", "Test")) %>%
  mutate(Task=factor(Task, levels=c("Training","Test","CB6F1", "B6-metagenomics"))) %>%
  ggplot(aes(x=Task, y=AUC)) +
  geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
  coord_cartesian(ylim=c(0,1), expand=0) +
  scale_fill_viridis_d()+theme_minimal()+ggtitle("AUC-Pathway")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18, color="black"), axis.text=element_text(size=18, color="black"),
        axis.text.x = element_text(angle = 30,hjust=1),axis.title.x = element_blank())


tiff("Fig_final/RF.ROC.pw.tiff", units = "in", height=6, width = 6.5, res=300, compression = "lzw")
roc.pw.b6
dev.off()

tiff("Fig_final/RF.AUC.pw.tiff", units = "in", height=6, width = 4.5, res=300, compression = "lzw")
auc.pw.b6
dev.off()

PW.logFC<-PW.full %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  as_tibble() %>%
  gather(-feature, key=id, value=Abundance) %>%
  left_join(RF.metadata %>% dplyr::select(id, Task, age)) %>%
  group_by(feature, age) %>%
  summarize(mean=mean(Abundance)) %>%
  spread(key=age, value=mean) %>%
  mutate(logFC=log2(old)-log2(young)) %>%
  ungroup() 

RF.features.pw.full<-RF.pw.full$rf.list[[10]]$importance %>% as.data.frame() %>% arrange(desc(MeanDecreaseGini)) %>% 
  top_n(20, MeanDecreaseGini) %>% rownames_to_column("feature") %>% left_join(PW.des, by='feature') %>% left_join(PW.logFC, by='feature')

tiff("Fig_final/RF.feature.pw.tiff", units = "in", height=8, width = 10, res=300, compression = "lzw")
RF.features.pw.full %>% arrange(desc(MeanDecreaseGini)) %>% 
  top_n(20, MeanDecreaseGini) %>% ggplot(aes(x=reorder(description,MeanDecreaseGini),y=MeanDecreaseGini,fill=logFC))+
  geom_bar(stat="identity") +theme_minimal() + coord_flip()+ scale_fill_gradient2(low="dodgerblue3", high="firebrick3",n.breaks=8) + 
  xlab("Pathway")+theme(text = element_text(size = 15), axis.text = element_text(color = "black"),axis.title.y = element_blank())+ggtitle("Age-specific Features")
dev.off()

write.csv(RF.features.pw.full,"PW.feature.csv")






















