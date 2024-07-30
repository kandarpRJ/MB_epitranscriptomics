library(WGCNA)
net<-readRDS("b.vsd_WGCNA_net.rds")
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

m6a_reg<-c("METTL3", "METTL14", "WTAP", "CBLL1", "VIRMA", "METTL16", "METTL5",
           "ZCCHC4", "TRMT112", "RBM15", "ZC3H13",
           "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "IGF2BP1",
           "IGF2BP2", "IGF2BP3", "FMR1",
           "FTO", "ALKBH5")

m6a_reg<-gl$ID[gl$SYMBOL%in%m6a_reg]
reg.cat<-c(rep("READER",3), "ERASER", "READER", rep("WRITER",3), "READER", "WRITER", "ERASER", rep("WRITER",2),
           rep("READER",2), rep("WRITER",5), rep("READER",2))

names(m6a_reg)<-reg.cat

m6aregmods<-data.frame(module=net$colors[names(net$colors)%in%m6a_reg])
m6aregmods$class<-unlist(lapply(rownames(m6aregmods), function(x) names(m6a_reg[m6a_reg==x])))
m6aregmods<-m6aregmods[order(m6aregmods$module),]
m6aregmods$modColor<-labels2colors(m6aregmods$module)

mod.sel<-names(net$colors[net$colors%in%m6aregmods$module])
mod.sel<-data.frame(genes=mod.sel) %>% tidyr::separate(genes, c("ID", "SYMBOL"))
genc.lnc<-read.table("gencode_v44.lncrna.txt",header = F, row.names = NULL)
mod.lnc.sel<-mod.sel[mod.sel$ID%in%genc.lnc$V1,]

library(foreach)
library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

ktest<-foreach(i=1:nrow(datlmplatform),.combine = rbind) %dopar% {
  l<-as.vector(t(datlmplatform[i,]))
  dat<-data.frame(exp=l,grp=b.meta$SUBGROUP)
  t<-rstatix::kruskal_test(dat, exp~grp)
  t$.y.<-rownames(datlmplatform)[i]
  t<-tidyr::separate(data = t, .y., sep = "_", c("gene_id", "gene_symbol"))
}
ktest$padj<-p.adjust(ktest$p, method = "bonferroni")

mod.lnc.sel<-merge(mod.lnc.sel, ktest, by.x="ID", by.y="gene_id")
mod.lnc.sel.sig<-mod.lnc.sel[mod.lnc.sel$padj<0.05,]
mod.lnc.sel.sig<-mod.lnc.sel.sig[order(mod.lnc.sel.sig$padj),]

library(survivalAnalysis)
library(survival)
library(survminer)
library(glmnet)
set.seed(10)
mcensor<-m[!is.na(m$OS_TIME)&!is.na(m$OS_STATUS)&m$SAMPLE%in%b.meta$SAMPLE,]
mcensor$OS_STATUS[mcensor$OS_TIME>10]<-0
mcensor$OS_TIME[mcensor$OS_TIME>10]<-10

x<-t(b.vsd[rownames(b.vsd)%in%paste(mod.lnc.sel$ID, mod.lnc.sel$gene_symbol, sep = "_"),mcensor$SAMPLE])
y<-Surv(mcensor$OS_TIME, mcensor$OS_STATUS)
fit<-glmnet(x, y, family="cox")
plot(fit)

library(caret)
set.seed(10)
flds <- createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(y))
foldids[flds$Fold02] = 2
foldids[flds$Fold03] = 3
foldids[flds$Fold04] = 4
foldids[flds$Fold05] = 5
foldids[flds$Fold06] = 6
foldids[flds$Fold07] = 7
foldids[flds$Fold08] = 8
foldids[flds$Fold09] = 9
foldids[flds$Fold10] = 10

cv.fit<-cv.glmnet(x, y, nfolds = 10, family="cox", foldid = foldids)
plot(cv.fit)

sel.coef<-as.matrix(coef(cv.fit, s="lambda.min"))
sel.coef<-sel.coef[sel.coef[,1]!=0,]
xy<-merge(x[,names(sel.coef)], y=b.meta[,c("SAMPLE", "OS_STATUS", "OS_TIME")], by.x=0, by.y="SAMPLE")
rownames(xy)<-xy$Row.names
xy<-xy[,-1]
xy$OS_STATUS[xy$OS_TIME>10]<-0
xy$OS_TIME[xy$OS_TIME>10]<-10

uni.cph.coef<-NULL
uni.cph.pval<-NULL
uni.cph.lowCI<-NULL
uni.cph.upCI<-NULL

for(i in 1:length(sel.coef)) {
  cph.res<-coxph(Surv(OS_TIME, OS_STATUS)~xy[,i], data=xy[,c(i,7,8)])
  names(cph.res$coefficients)<-names(sel.coef)[i]
  uni.cph.coef<-c(uni.cph.coef, cph.res$coefficients)
  uni.cph.pval<-c(uni.cph.pval, summary(cph.res)$coefficient[5])
  uni.cph.lowCI<-c(uni.cph.lowCI, summary(cph.res)$conf.int[3])
  uni.cph.upCI<-c(uni.cph.upCI, summary(cph.res)$conf.int[4])
}

rem<-which(uni.cph.pval>=0.05)
uni.cph.coef<-uni.cph.coef[-rem]
uni.cph.pval<-uni.cph.pval[-rem]
uni.cph.lowCI<-uni.cph.lowCI[-rem]
uni.cph.upCI<-uni.cph.upCI[-rem]
uni.cph<-data.frame(coef=uni.cph.coef, pval=uni.cph.pval, HR=exp(uni.cph.coef), Lower_CI=uni.cph.lowCI, Upper_CI=uni.cph.upCI)
uni.cph<-uni.cph[order(uni.cph$HR, decreasing = T),]

ggplot(uni.cph, aes(HR, reorder(gsub("ENSG\\w*_", "",rownames(uni.cph)),HR), color = HR > 1)) +
  geom_vline(xintercept = 1, color = "gray75") +
  geom_linerange(aes(xmin = Lower_CI, xmax = Upper_CI), size = 1.5, alpha = 0.5) +
  geom_point(size = 4) +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = c("green4", "red3"), guide = "none") +
  xlim(c(0, 2.5)) +
  labs(title = "Hazard ratio for m6A-lncRNA signature genes", y = NULL,
       x = "Hazard ratio estimate (95% Confidence Intervals)") +
  theme(axis.text.y = element_text(hjust = 0, size = 18))

uni.cph<-data.frame(coef=uni.cph.coef, pval=uni.cph.pval)

library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
cr<-read.csv("cellrep.raw.counts.mb.csv")
cr$X<-gsub("\\.\\w*", "", cr$X)

gl<-clusterProfiler::bitr(cr$X, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop=F)
for (i in 1:nrow(gl)) {
  if (is.na(gl$SYMBOL[i])) {
    gl$SYMBOL[i]=gl$ENSEMBL[i]
  }
}
gl$ID<-paste(gl$ENSEMBL, gl$SYMBOL, sep="_")

cr<-merge(gl, cr, by.x="ENSEMBL", by.y="X")
cr<-cr %>%
  mutate(Gene_Symbol = paste(ENSEMBL, SYMBOL, sep="_"))
cr<-cr[,c("Gene_Symbol",colnames(cr)[4:ncol(cr)-1])]
cr<-cr%>%
  group_by(Gene_Symbol) %>%
  summarise(across(where(is.numeric), max))
cr<-data.frame(cr)
rownames(cr)<-cr$Gene_Symbol
cr<-cr[,-1]

crmeta<-read.csv(gzfile("cellrep_meta.csv.gz"))
crmeta$SUBGROUP[crmeta$SUBGROUP=="G3"]<-"Grp3"
crmeta$SUBGROUP[crmeta$SUBGROUP=="G4"]<-"Grp4"

pbta<-read.csv(gzfile("pbta_counts.csv.gz"), row.names = 1)
pbtameta<-read.csv(gzfile("pbta_meta.csv.gz"))
pbtameta$SUBGROUP[pbtameta$SUBGROUP=="G3"]<-"Grp3"
pbtameta$SUBGROUP[pbtameta$SUBGROUP=="G4"]<-"Grp4"

stjude<-read.csv(gzfile("stjude_counts.csv.gz"), row.names = 1)
stjudemeta<-read.csv(gzfile("stjude_meta.csv.gz"))
stjudemeta$SUBGROUP[stjudemeta$SUBGROUP=="G3"]<-"Grp3"
stjudemeta$SUBGROUP[stjudemeta$SUBGROUP=="G4"]<-"Grp4"

library(GEOquery)
gse<-getGEO("GSE85217", GSEMatrix = T, AnnotGPL = T)
exp<-exprs(gse[[1]])
exp<-tibble::rownames_to_column(as.data.frame(exp))
exp$rowname<-gsub("_at", "", exp$rowname)
exp<-merge(gl, exp, by.x="ENSEMBL", by.y="rowname")
exp<-exp[,c(3:ncol(exp))]

gse.meta<-read.csv(gzfile("gse_meta.csv.gz"))
gse.meta$SUBGROUP[gse.meta$SUBGROUP=="G3"]<-"Grp3"
gse.meta$SUBGROUP[gse.meta$SUBGROUP=="G4"]<-"Grp4"

rownames(exp)<-exp$ID
exp<-exp[,-1]

m<-rbind(crmeta, pbtameta, stjudemeta, gse.meta)
m<-m[m$SUBGROUP!="MBL",]
m<-m[m$SUBGROUP!="Unknown",]
plat<-data.frame(sample_name=m$SAMPLE, Class=c(rep("1", 308), rep("2", 102), rep("3",89), rep("4", 737)))
m$PLATFORM<-plat$Class
m<-m[!is.na(m$AGE)&!is.na(m$SEX),]
m$AGE[m$AGE=="0-3"]<-1
m$AGE[m$AGE=="3-16"]<-2
m$AGE[m$AGE=="over16"]<-3
m$SOUT<-1
m$SOUT[m$SUBGROUP=="SHH"]<-2
m$SOUT[m$SUBGROUP=="Grp3"]<-3
m$SOUT[m$SUBGROUP=="Grp4"]<-4

library(DESeq2)

a<-merge(cr, pbta, by=0)
b<-merge(a, stjude, by.x="Row.names", by.y=0)
b<-b[b$Row.names%in%rownames(exp),]
rownames(b)<-b$Row.names
b<-b[,-1]
msam<-m$SAMPLE[m$SAMPLE%in%colnames(b)]
b.meta<-m[m$SAMPLE%in%msam,]
b<-b[,msam]
dds<-DESeqDataSetFromMatrix(b, m[m$SAMPLE%in%msam,], ~SUBGROUP)
dds<-dds[rowSums(counts(dds)>10)>234,]
dds<-DESeq(dds)
vsd<-varianceStabilizingTransformation(dds, blind = T)
plotPCA(vsd, intgroup="SUBGROUP")+theme_light()+theme(text=element_text(size = 14))
b.vsd<-assay(vsd)
b.vsd<-round(b.vsd,2)
sam.merge<-b.vsd
risk_score<-data.frame(apply(sam.merge[rownames(sam.merge)%in%rownames(uni.cph),]*uni.cph$coef,2,sum))
colnames(risk_score)<-"risk_score"

mcensor<-m[!is.na(m$OS_STATUS)&!is.na(m$OS_TIME),]
mcensor$OS_STATUS[mcensor$OS_TIME>10]<-0
mcensor$OS_TIME[mcensor$OS_TIME>10]<-10

risk_score<-merge(risk_score, mcensor, by.x=0, by.y="SAMPLE")
risk_score$AGE<-as.numeric(risk_score$AGE)
risk_score$OS_TIME<-risk_score$OS_TIME*12

risk_score$SOUT<-1
risk_score$SOUT[risk_score$SUBGROUP=="SHH"]<-2
risk_score$SOUT[risk_score$SUBGROUP=="G3"]<-3
risk_score$SOUT[risk_score$SUBGROUP=="G4"]<-4

risk_score<-risk_score[!is.na(risk_score$OS_TIME)&!is.na(risk_score$OS_STATUS),]
risk_score$group<-kmeans(risk_score$risk_score, centers = 2)$cluster

med.exp<-median(risk_score$risk_score)
more.med.exp.index<-which(risk_score$risk_score>=med.exp)
less.med.exp.index<-which(risk_score$risk_score< med.exp)
risk_score$status<-NA
risk_score$status[more.med.exp.index]<-'High'
risk_score$status[less.med.exp.index]<-'Low'

sub.risk_score<-risk_score[risk_score$PLATFORM==2,]
sub.risk_score$group<-kmeans(sub.risk_score$risk_score, centers = 2)$cluster

med.exp<-median(sub.risk_score$risk_score)
more.med.exp.index<-which(sub.risk_score$risk_score>=med.exp)
less.med.exp.index<-which(sub.risk_score$risk_score< med.exp)
sub.risk_score$status<-NA
sub.risk_score$status[more.med.exp.index]<-'High'
sub.risk_score$status[less.med.exp.index]<-'Low'

res<-coxph(Surv(OS_TIME, OS_STATUS) ~ status + AGE + SEX, data=risk_score)
ggforest(res, fontsize = 1.2)

s.fit<-survfit(Surv(OS_TIME, OS_STATUS) ~ status, data = risk_score)
s.diff<-survdiff(Surv(OS_TIME, OS_STATUS) ~ status, data = risk_score)

sdata.plot1<-ggsurvplot(s.fit,
                        data=risk_score,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (years)',
                        ggtheme = ggplot2::theme_light(base_size = 14),
                        surv.median.line = 'hv',
                        title="OS [combined]")

sdata.plot1

s.fit<-survfit(Surv(OS_TIME, OS_STATUS) ~ status, data = sub.risk_score)
s.diff<-survdiff(Surv(OS_TIME, OS_STATUS) ~ status, data = sub.risk_score)

sdata.plot1<-ggsurvplot(s.fit,
                        data=sub.risk_score,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (years)',
                        ggtheme = ggplot2::theme_light(base_size = 14),
                        surv.median.line = 'hv',
                        title="OS [Study 2]")
sdata.plot1

library(rms)
set.seed(10)
num <- 10
fold=1
indices <- sample(1:num, size = nrow(risk_score), replace = TRUE)

train_data <- risk_score[indices != fold,]
test_data <- risk_score[indices == fold,]

ddist <- datadist(train_data)
options(datadist='ddist')

f <- cph(Surv(OS_TIME, OS_STATUS) ~ risk_score + AGE + SEX,
         x=T, y=T, surv=T, data=train_data, time.inc=3)

surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1*12, x),
                            function(x) surv(3*12, x), function(x) surv(5*12, x)), lp=F, 
                funlabel=c("1-year survival", "3-year survival", "5-year survival"),
                risk_score=seq(-13,-7, by=0.1), AGE=c(1,2,3), SEX=c(1,2), fun.at=seq(1,0,by=-0.1))

plot(nom, cex.main=2.5, cex.lab=2.5)

nom <- nomogram(f, fun=list(function(x) surv(1*12, x),
                            function(x) surv(3*12, x), function(x) surv(5*12, x)), lp=T, 
                funlabel=c("1-year survival", "3-year survival", "5-year survival"),
                risk_score=seq(-13,-7, by=0.1), AGE=c(1,2,3), SEX=c(1,2), fun.at=seq(1,0,by=-0.1))

library(pROC)
lps<-predict(f, test_data, type = "lp")
predicted_probs<-plogis(lps)
risk_groups<- cut(predicted_probs, quantile(predicted_probs, probs = seq(0, 1, 0.25)), 
                  include.lowest = TRUE, labels = FALSE)
roc_curve <- roc(test_data$OS_STATUS, predicted_probs)
plot(roc_curve, main = "ROC Curve", col = "black")
text(x=0,y=0.1, labels=paste("AUC: ", round(roc_curve$auc[1],2)), cex=1.2, pos=2)

f1 <- cph(Surv(OS_TIME, OS_STATUS) ~ risk_score + AGE + SEX,
          x=T, y=T, surv=T, data=risk_score, time.inc=1)
cal<-rms::calibrate(f1,cmethod="KM", method="boot", u=1, xlim=c(0,1), ylim=c(0,1), m=50)
plot(cal, xlim=c(0,1), ylim=c(0,1), col="green4")

f3 <- cph(Surv(OS_TIME, OS_STATUS) ~ risk_score + AGE + SEX,
          x=T, y=T, surv=T, data=risk_score, time.inc=3)
cal<-rms::calibrate(f3,cmethod="KM", method="boot", u=3, xlim=c(0,1), ylim=c(0,1), m=50)
plot(cal, xlim=c(0,1), ylim=c(0,1), add=T, col="blue4")

f5 <- cph(Surv(OS_TIME, OS_STATUS) ~ risk_score + AGE + SEX,
          x=T, y=T, surv=T, data=risk_score, time.inc=5)
cal<-rms::calibrate(f5,cmethod="KM", method="boot", u=5, xlim=c(0,1), ylim=c(0,1), m=50, add=T)
plot(cal, xlim=c(0,1), ylim=c(0,1), add=T, col="red4")

legend('bottomright',
       legend = c('1-year survival',
                  '3-year survival',
                  '5-year survival'),
       col=c('blue4','red4', 'green4'),lwd=2, cex = 2)

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

asian.c<-read.table("asian_mb_counts.tsv", sep = "\t", header = T, row.names = 1)
x<-clusterProfiler::bitr(rownames(asian.c), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
x$id<-paste(x$ENSEMBL, x$SYMBOL, sep="_")
asian.c<-merge(asian.c,x,by.x=0,by.y="ENSEMBL")
rownames(asian.c)<-asian.c$id
asian.c<-asian.c[,-c(1,61,62)]

asian_meta<-read.table("asian_mb_meta.tsv", sep = "\t", header = T)
rownames(asian_meta)<-asian_meta$RUN
asian_meta<-asian_meta[colnames(asian.c),]
asian_meta$SOUT<-1
asian_meta$SOUT[asian_meta$SUBGROUP=="MB-SHH"]<-2
asian_meta$SOUT[asian_meta$SUBGROUP=="MB-G3"]<-3
asian_meta$SOUT[asian_meta$SUBGROUP=="MB-G4"]<-4
asian_meta$PLATFORM<-5

brazil.c<-read.table("brazil_mb_counts.tsv", sep="\t", header = T, row.names = 1)
x<-clusterProfiler::bitr(rownames(brazil.c), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
x$id<-paste(x$ENSEMBL, x$SYMBOL, sep="_")
brazil.c<-merge(brazil.c,x,by.x=0,by.y="ENSEMBL")
rownames(brazil.c)<-brazil.c$id
brazil.c<-brazil.c[,-c(1,19,20)]
bm<-read.table("brazil_mb_meta.tsv", header = T, sep = "\t")
bm$PLATFORM<-6

sb.trin<-merge(b, asian.c, by=0)
sb.trin<-merge(sb.trin, brazil.c, by.x="Row.names", by.y=0)
rownames(sb.trin)<-sb.trin$Row.names
sb.trin<-sb.trin[,-1]

dds<-DESeqDataSetFromMatrix(sb.trin, 
                            colData = data.frame(SOUT=c(b.meta$SOUT, asian_meta$SOUT, bm$SOUT),
                                                 STUDY=c(b.meta$PLATFORM, asian_meta$PLATFORM, bm$PLATFORM)), ~SOUT+STUDY)
dds<-DESeq(dds)
vsd<-assay(varianceStabilizingTransformation(dds, blind = F))

asian_meta$SUBGROUP2<-"nonWNT_nonSHH"
asian_meta$SUBGROUP2[asian_meta$SUBGROUP=="WNT"]<-"WNT"
asian_meta$SUBGROUP2[asian_meta$SUBGROUP=="SHH"]<-"SHH"
bm$SUBGROUP2<-"nonWNT_nonSHH"
bm$SUBGROUP2[bm$SUBGROUP=="WNT"]<-"WNT"
bm$SUBGROUP2[bm$SUBGROUP=="SHH"]<-"SHH"
meta.df<-data.frame(samples=c(b.meta$SAMPLE, asian_meta$SAMPLE_NAME, bm$run_accession), subgroups=c(b.meta$SUBGROUP2, asian_meta$SUBGROUP2, bm$SUBGROUP2), 
                    STUDY=c(b.meta$PLATFORM, asian_meta$PLATFORM, bm$PLATFORM))

m6a_reg<-m6a_reg[order(names(m6a_reg))]
library(ComplexHeatmap)
hm.mat<-b.vsd[rownames(b.vsd)%in%paste(mod.lnc.sel$ID, mod.lnc.sel$gene_symbol, sep = "_"),]
hm.mat<-t(scale(t(hm.mat)))
rownames(hm.mat)<-gsub("\\w*_","",rownames(hm.mat))
ha<-HeatmapAnnotation(Subgroups=b.meta[,c("SUBGROUP2")], col = list("Subgroups"=c("nonWNT_nonSHH"="red", "SHH"="green", "WNT"="black"),
                                                                                      "PLATFORM"=c("1"="gold", "2"="purple", "3"="cyan", "5"="pink", "6"="magenta")))
dend<-cluster_within_group(hm.mat,b.meta$SUBGROUP2)
Heatmap(hm.mat, show_column_names = F, cluster_rows = T, cluster_columns = dend,
        top_annotation = ha, name = "Scaled expression")

all_merged<-merge(vsd, exp, by=0)
rownames(all_merged)<-all_merged$Row.names
all_merged<-all_merged[,-1]


gse.meta$SOUT<-1
gse.meta$SOUT[gse.meta$SUBGROUP=="SHH"]<-2
gse.meta$SOUT[gse.meta$SUBGROUP=="G3"]<-3
gse.meta$SOUT[gse.meta$SUBGROUP=="G4"]<-4

gse.meta$PLATFORM<-4

all.meta<-rbind(b.meta[,c("SOUT", "PLATFORM")], asian_meta[,c("SOUT", "PLATFORM")],
                bm[,c("SOUT", "PLATFORM")], gse.meta[,c("SOUT", "PLATFORM")])
rownames(all.meta)<-c(b.meta$SAMPLE, asian_meta$RUN, bm$run_accession, gse.meta$SAMPLE)
all.meta<-all.meta[colnames(all_merged),]

com_corr<-sva::ComBat(dat=all_merged, batch = all.meta$PLATFORM)

all.pr<-prcomp(t(com_corr), scale. = T)
ggplot(as.data.frame(all.pr$x), ggplot2::aes(PC1, PC2, col=as.factor(all.meta$SOUT)))+ggplot2::geom_point()

sam.merge<-com_corr[,b.meta$SAMPLE]

bor.exp<-t(sam.merge[c(paste(mod.lnc.sel$ID, mod.lnc.sel$gene_symbol, sep = "_")),])
bor.grp<-data.frame(subgroup=b.meta$SOUT)
rownames(bor.grp)<-b.meta$SAMPLE
train.dat<-merge(bor.exp,bor.grp,by=0)
train.dat<-na.omit(train.dat)

library(Boruta)
set.seed(10)
boruta.train <- Boruta(factor(subgroup)~.-Row.names-subgroup, data = train.dat, doTrace = 2, maxRuns=100)
print (boruta.train)
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
fb<-final.boruta
colnames(fb$ImpHistory)<-gsub("\\w*_","",colnames(fb$ImpHistory))
par(mar=c(6,2,2,2))
plot(fb, las=2, cex.axis=0.8, xlab="")

trainin<-com_corr[,b.meta$SAMPLE]
trainin<-trainin[names(final.boruta$finalDecision[final.boruta$finalDecision=="Confirmed"]),b.meta$SAMPLE]
trainin<-cbind(t(trainin), data.frame(SOUT=b.meta$SOUT))
train_id<-createDataPartition(trainin$SOUT, p=0.8, list = F)
train_data<-trainin[train_id,]
test_data<-trainin[-train_id,]

trainControl <- trainControl(method="repeatedcv", number=5, repeats=3, verboseIter = TRUE)

metric <- "Accuracy"
set.seed(10)
fit.rf <- train(as.factor(SOUT)~., data = train_data, method = "rf", metric = metric,trControl = trainControl)
fit.gbm <- train(as.factor(SOUT)~., data = train_data, method = "gbm",metric = metric,trControl = trainControl, verbose = FALSE)
fit.c50 <- train(as.factor(SOUT)~., data = train_data, method = "C5.0", metric = metric,trControl = trainControl)
fit.lda <- train(as.factor(SOUT)~., data = train_data, method="lda",
                 metric=metric,trControl=trainControl)
fit.glmnet <- train(as.factor(SOUT)~., data = as.matrix(train_data), method="glmnet",
                    metric=metric,trControl=trainControl)
fit.knn <- train(as.factor(SOUT)~., data = train_data, method="knn",
                 metric=metric,trControl=trainControl)
fit.svm <- train(as.factor(SOUT)~., data = train_data, method="svmRadial",
                 metric=metric,trControl=trainControl)
fit.xgboost<-train(as.factor(SOUT)~., data = train_data, method="xgbTree", 
                   metric=metric, trControl=trainControl)
results <- resamples(list(RF=fit.rf, GBM=fit.gbm, C50=fit.c50, LDA=fit.lda, 
                          GLMNET=fit.glmnet, KNN=fit.knn,  
                          SVM=fit.svm)) #, XGBOOST=fit.xgboost))
summary(results)
bwplot(results)
dotplot(results)

pred<-predict(fit.glmnet, test_data[,1:67],probability=T)
confusionMatrix(pred, as.factor(test_data$SOUT), mode="everything")
library(pROC)
mroc<-pROC::multiclass.roc(as.numeric(test_data$SOUT)~as.numeric(pred), plot=T, print.auc=T)
plot.roc(mroc$rocs[[1]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5)
plot.roc(mroc$rocs[[2]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5, add = T, col = "blue", print.auc.adj = c(0,3))
plot.roc(mroc$rocs[[3]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5, add = T, col = "red", print.auc.adj = c(0,5))
plot.roc(mroc$rocs[[4]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5, add = T, col = "darkgreen", print.auc.adj = c(0,7))
plot.roc(mroc$rocs[[5]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5, add = T, col = "orange", print.auc.adj = c(0,9))
plot.roc(mroc$rocs[[6]], print.auc = T, legacy.axes = T, print.auc.cex = 1.5, add = T, col = "purple", print.auc.adj = c(0,11))
legend('bottomright',
       legend = c('WNT-SHH',
                  'WNT-G3',
                  'WNT-G4',
                  'SHH-G3',
                  'SHH-G4',
                  'G3-G4'),
       col=c('black','blue','red', 'darkgreen', 'orange', 'purple'),lwd=2, cex = 2)

celltype<-read.csv("CIBERSORTx_Job12_Results.csv")
celltype<-celltype[celltype$P.value<0.05,]

samples<-celltype$Mixture[celltype$Mixture%in%colnames(sam.merge)]

cel.cor<-foreach(i = 1:nrow(uni.cph),.combine=rbind) %dopar% {
  gexp<-as.vector(t(sam.merge[rownames(sam.merge)%in%rownames(uni.cph)[i],samples]))
  x=data.frame()
  for (j in 2:23) {
    tres<-cor.test(celltype[celltype$Mixture%in%samples,j],gexp)
    x<-rbind(x, data.frame(colnames(celltype)[j],rownames(uni.cph)[i],tres$estimate, tres$p.value))
  }
  x
}
cel.cor.sig<-cel.cor[cel.cor$tres.p.value<0.05,]
cel.cor.sig<-cel.cor.sig[order(cel.cor.sig$tres.estimate),]
colnames(cel.cor.sig)<-c("Cell_type", "lncRNA", "cor", "pvalue")
cel.cor.sig$lncRNA<-gsub("\\w*_","",cel.cor.sig$lncRNA)

par(mar = c(2, 2, 2, 2))

ggplot2::ggplot(cel.cor.sig, aes(Cell_type, lncRNA, col=cor, size=-log(pvalue))) +
  geom_point() +
  xlab("Cell types")+
  scale_color_gradient(low="darkblue", high = "gold") +
  theme_light() +
  theme(plot.margin=unit(c(1,1,2,1), units = "cm"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

res=NULL
sub.celltype<-celltype[celltype$Mixture%in%risk_score$Row.names,]
for (k in 2:23) {
  ct.cor<-cor.test(risk_score$risk_score[risk_score$Row.names%in%sub.celltype$Mixture], sub.celltype[,k])
  res<-rbind(res, data.frame(colnames(celltype)[k],ct.cor$estimate, ct.cor$p.value))
}

b.meta.rs<-risk_score
rownames(b.meta.rs)<-b.meta.rs$Row.names
b.meta.rs<-b.meta.rs[,-1]
b.meta.rs$status<-as.factor(b.meta.rs$status)
b.meta.rs$status<-relevel(b.meta.rs$status, ref = "Low")
b.rs<-b[,rownames(b.meta.rs)]
dds<-DESeqDataSetFromMatrix(b.rs, b.meta.rs, ~status)
dds<-DESeq(dds)
deg.res<-results(dds, alpha = 0.05)
deg.res<-data.frame(deg.res)
deg.res<-deg.res[deg.res$padj<0.05,]
deg.res<-deg.res[order(deg.res$padj, decreasing = F),]

volcano_plot <- function(res, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  sig_res <- subset(res, padj <= padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(color = ifelse(abs(res$log2FoldChange) >= log2fc_cutoff & res$padj <= padj_cutoff, "red", "black"), alpha = 0.6) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "blue") +
    labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)", title = "Volcano Plot") +
    theme_minimal() +
    theme(legend.position="none")
}

volcano_plot(deg.res)

library(clusterProfiler)
library(enrichplot)
geneli<-deg.res$stat
names(geneli)<-gsub("\\w*_", "", rownames(deg.res))
geneli<-sort(geneli, decreasing = T)

gobp<-enrichGO(gene = names(geneli), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
dotplot(gobp)

ego<-gseGO(geneList = geneli, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
gseaplot2(ego, geneSetID = 10, title = ego$Description[10])

icgenes<-c("CD276", "CTLA4", "CD80", "CD86", "PDCD1", "CD274", "PDCD1LG2", "LAG3", "CD158", "HAVCR2", "CEACAM1", 
           "HMGB1", "TIGIT", "CD226", "CD96", "PVR",  "TNFRSF9", "TNFSF9", "TNFRSF18", "TNFSF18")

icgenes.sel<-rownames(b.rs)[unlist(lapply(icgenes, function (x) grep (x, rownames(b.rs))))]
icgenes.sel<-icgenes.sel[c(1:4,7,9:13,16:24)]

names(icgenes.sel)<-c("Ligand", "Receptor", "Ligand", "Ligand", "Receptor", "Ligand", "Ligand", "Receptor", "Ligand", "Ligand", "Receptor", "Receptor", "Receptor", "Ligand", "Receptor", "Ligand", "Receptor", "Ligand", "Receptor")

de.icgenes<-na.omit(deg.res[icgenes.sel,])

b.rs.vsd<-assay(varianceStabilizingTransformation(dds)[rownames(dds)%in%icgenes.sel,])
icgenes.label<-rownames(b.rs.vsd)
icgenes.label<-gsub("\\w*_","",icgenes.label)
icgenes.label<-c("4-1BB","PVR","CEACAM1","LAG3","CD276","CD86","CD274","HGITRL",
                 "CD80","4-1BBL","HAVCR2","CD226","CD96","CTLA4","TIGIT","GITR","CD279","HMGB1","CD273")
names(icgenes.label)<-c(rep("Receptor", 4), rep("Ligand", 6), rep("Receptor", 7), rep("Ligand", 2))
row.ha<-rowAnnotation(Gene_type=names(icgenes.label), 
                      col=list(Gene_type=c("Ligand"="gold","Receptor"="steelblue")))
ha.wnt<-HeatmapAnnotation(subgroup=b.meta.rs$SUBGROUP[b.meta.rs$SUBGROUP=="WNT"], 
                          risk_score=b.meta.rs$status[b.meta.rs$SUBGROUP=="WNT"], show_annotation_name = F,
                          col = list(subgroup=c("WNT"="black"), risk_score=c("Low"="cyan", "High"="magenta")))
wnt.mat<-as.matrix(t(scale(t(b.rs.vsd[,rownames(b.meta.rs)[b.meta.rs$SUBGROUP=="WNT"]]))))
rownames(wnt.mat)<-gsub("\\w*_","",rownames(wnt.mat))
p.wnt<-Heatmap(wnt.mat, show_column_names = F, cluster_rows = F,
               top_annotation = ha.wnt, name = "Z-score")

ha.shh<-HeatmapAnnotation(subgroup=b.meta.rs$SUBGROUP[b.meta.rs$SUBGROUP=="SHH"], 
                          risk_score=b.meta.rs$status[b.meta.rs$SUBGROUP=="SHH"], show_annotation_name = F,
                          col = list(subgroup=c("SHH"="green"), risk_score=c("Low"="cyan", "High"="magenta")))
shh.mat<-as.matrix(t(scale(t(b.rs.vsd[,rownames(b.meta.rs)[b.meta.rs$SUBGROUP=="SHH"]]))))
rownames(shh.mat)<-gsub("\\w*_","",rownames(shh.mat))
p.shh<-Heatmap(shh.mat, show_column_names = F, cluster_rows = F, 
               top_annotation = ha.shh, name = "Z-score")

ha.grp3<-HeatmapAnnotation(subgroup=b.meta.rs$SUBGROUP[b.meta.rs$SUBGROUP=="Grp3"], 
                           risk_score=b.meta.rs$status[b.meta.rs$SUBGROUP=="Grp3"], show_annotation_name = F,
                           col = list(subgroup=c("Grp3"="red"), risk_score=c("Low"="cyan", "High"="magenta")))
grp3.mat<-as.matrix(t(scale(t(b.rs.vsd[,rownames(b.meta.rs)[b.meta.rs$SUBGROUP=="Grp3"]]))))
rownames(grp3.mat)<-gsub("\\w*_","",rownames(grp3.mat))
p.grp3<-Heatmap(grp3.mat, show_column_names = F, cluster_rows = F, 
                top_annotation = ha.grp3, name = "Z-score")

ha.grp4<-HeatmapAnnotation(subgroup=b.meta.rs$SUBGROUP[b.meta.rs$SUBGROUP=="Grp4"], 
                           risk_score=b.meta.rs$status[b.meta.rs$SUBGROUP=="Grp4"], 
                           col = list(subgroup=c("Grp4"="blue"), risk_score=c("Low"="cyan", "High"="magenta")))
grp4.mat<-as.matrix(t(scale(t(b.rs.vsd[,rownames(b.meta.rs)[b.meta.rs$SUBGROUP=="Grp4"]]))))
rownames(grp4.mat)<-gsub("\\w*_","",rownames(grp4.mat))
p.grp4<-Heatmap(grp4.mat, show_column_names = F, cluster_rows = F, 
                right_annotation = row.ha, row_labels = icgenes.label, top_annotation = ha.grp4, name = "Z-score")

p.wnt+p.shh+p.grp3+p.grp4

d1<-data.frame(gene="PVR", gene_exp=b.rs.vsd["ENSG00000073008_PVR",], status=b.meta.rs$status)
colnames(d1)<-c("gene", "gene_exp", "status")
d2<-data.frame(gene="4-1BB", gene_exp=b.rs.vsd["ENSG00000049249_TNFRSF9",], status=b.meta.rs$status)
colnames(d2)<-c("gene", "gene_exp", "status")
d3<-data.frame(gene="4-1BBL", gene_exp=b.rs.vsd["ENSG00000125657_TNFSF9",], status=b.meta.rs$status)
colnames(d3)<-c("gene", "gene_exp", "status")
d4<-data.frame(gene="GITR", gene_exp=b.rs.vsd["ENSG00000186891_TNFRSF18",], status=b.meta.rs$status)
colnames(d4)<-c("gene", "gene_exp", "status")
d5<-data.frame(gene="HMGB1", gene_exp=b.rs.vsd["ENSG00000189403_HMGB1",], status=b.meta.rs$status)
colnames(d5)<-c("gene", "gene_exp", "status")
d6<-data.frame(gene="CD273", gene_exp=b.rs.vsd["ENSG00000197646_PDCD1LG2",], status=b.meta.rs$status)
colnames(d6)<-c("gene", "gene_exp", "status")
d<-rbind(d1, d2, d3, d4, d5, d6)

ggplot(d, aes(gene, gene_exp))+
  geom_boxplot(aes(fill=status)) +
  facet_wrap(~gene, scales = "free") +
  xlab("") +
  ylab("Normalized expression") +
  theme_light() +
  theme(text=element_text(size=24))

