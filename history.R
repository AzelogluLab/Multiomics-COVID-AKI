

##############################
#feature selection 
library(Boruta)
set.seed(111)


genes<-read.delim(paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p.txt", sep=''), check.names = F)
sample_info<-read.delim("covid_proteomics_samples_sum.txt")
sample_info<-sample_info[sample_info$Set=="train"&!is.na(sample_info$Composite_Outcome),]
sample_info$Study_ID<-sample_info$Sample

genes_mat<-genes[, as.character(sample_info$Study_ID)]
rownames(genes_mat)<-genes$Accession
genes_mat<-t(genes_mat)

genes_mat<-merge(genes_mat, sample_info[,c("Study_ID", "Composite_Outcome")], by.x="row.names", by.y="Study_ID")
rownames(genes_mat)=genes_mat[,1]
genes_mat<-genes_mat[,-1]



feature_sel_train <- Boruta(Composite_Outcome~., data = genes_mat, doTrace = 2, maxRuns=1000)
feature_sel <- TentativeRoughFix(feature_sel_train)

print(feature_sel)

index<-which(feature_sel$finalDecision=="Confirmed")


genes<-read.delim(paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p.txt", sep=''), check.names = F)
#annot<-genes[, 1:6]
genes$Symbol<-gsub(".*GN=", "", genes$Description)
genes$Symbol<-gsub(" .*", "", genes$Symbol)

lz<-lapply(1:ncol(feature_sel$ImpHistory),function(i)feature_sel$ImpHistory[is.finite(feature_sel$ImpHistory[,i]),i])
names(lz) <- colnames(feature_sel$ImpHistory)
lz_cofirmed<-lz[names(feature_sel$finalDecision)[index]]

Labels <- sort(sapply(lz_cofirmed,median))
lz_cofirmed<-lz_cofirmed[names(Labels)]

gn<-merge(data.frame(Accession=names(lz_cofirmed)),genes[,c(1,2)], by="Accession", sort=F)

lz_cofirmed_df <- data.frame(Importance = unlist(lz_cofirmed), 
                             Gene = rep(as.character(gn$Symbol),times = sapply(lz_cofirmed,length)))
lz_cofirmed_df$Gene<-factor(lz_cofirmed_df$Gene, levels=as.character(gn$Symbol))

library(ggplot2)
pdf(paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_selected_features_importance.pdf", sep=""), width=4, height=5)
ggplot(lz_cofirmed_df, aes(x=Gene, y=Importance)) +
  geom_boxplot(color="black")+coord_flip()+theme_classic()
dev.off()

write.table(lz_cofirmed_df,  paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_selected_features_importance.txt", sep=''), row.names=F, sep="\t", na="", quote=F)


genes_sel<-merge(data.frame(Accession=names(Labels), med_imp=Labels), genes, by="Accession", sort=F)
write.table(genes_sel,  paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel.txt", sep=''), row.names=F, sep="\t", na="", quote=F)



#####################
#random forest 10-f CV adj
library(randomForest)
set.seed(70)
library(caret)
#read dataset
genes_sel<-read.delim(paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel.txt", sep=''), check.names = F)
genes_sel<-genes_sel[order(genes_sel$med_imp, decreasing = T), ]
sample_info<-read.delim("covid_proteomics_samples_sum.txt")
sample_info<-sample_info[sample_info$Set=="train"&!is.na(sample_info$Composite_Outcome),]
sample_info$Study_ID<-sample_info$Sample


genes_mat<-genes_sel[, as.character(sample_info$Study_ID)]
rownames(genes_mat)<-genes_sel$Accession
genes_mat<-t(genes_mat)

genes_mat<-merge(genes_mat, sample_info[,c("Study_ID", "Composite_Outcome")], by.x="row.names", by.y="Study_ID")
rownames(genes_mat)=genes_mat[,1]
genes_mat<-genes_mat[,-1]
#genes_mat$Composite_Outcome<-as.factor(genes_mat$Composite_Outcome)
genes_mat$Class<-rep(NA, nrow(genes_mat))
genes_mat$Class[genes_mat$Composite_Outcome==1]="Y"
genes_mat$Class[genes_mat$Composite_Outcome==0]="N"
genes_mat$Class<-as.factor(genes_mat$Class)
genes_mat<-genes_mat[, -(ncol(genes_mat)-1)]



#CV
#tune mtry
set.seed(70)
#create tunegrid
tunegrid <- expand.grid(.mtry = 1:(ncol(genes_mat)-1))
rf_gridsearch <- train(Class ~ ., 
                       data = genes_mat,
                       method = 'rf',
                       metric = 'Accuracy',
                       tuneGrid = tunegrid)
print(rf_gridsearch)
bestmtry<-as.numeric(rf_gridsearch$bestTune)


#tune ntree
control <- trainControl(method = 'repeatedcv',
                        number = 10,
                        repeats = 3,
                        search = 'grid',classProbs = TRUE)
#create tunegrid
tunegrid <- expand.grid(.mtry = bestmtry)
modellist <- list()

for (ntree in c(500, 1000,1500,2000,2500)){
  set.seed(70)
  fit <- train(Class~.,
               data = genes_mat,
               method = 'rf',
               metric = 'Accuracy',
               tuneGrid = tunegrid,
               trControl = control,
               ntree = ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}
results <- resamples(modellist)
a<-summary(results)
best_ntree<-as.numeric(names(which.max(a$statistics$Accuracy[,4])))

#retune mtry
if(best_ntree!=500){
  #create tunegrid
  tunegrid <- expand.grid(.mtry = 1:(ncol(genes_mat)-1))
  rf_gridsearch <- train(Class ~ ., 
                         data = genes_mat,
                         method = 'rf',
                         metric = 'Accuracy',
                         tuneGrid = tunegrid,
                         ntree = best_ntree)
  print(rf_gridsearch)
  bestmtry<-as.numeric(rf_gridsearch$bestTune)
}


#best model
set.seed(70)
control <- trainControl(method = 'repeatedcv',
                        number = 10,
                        repeats = 3,
                        search = 'grid',classProbs = TRUE,savePredictions="final")
tunegrid <- expand.grid(.mtry = bestmtry)
best_fit<-train(Class~.,
                data = genes_mat,
                method = 'rf',
                metric = 'Accuracy',
                tuneGrid = tunegrid,
                trControl = control,
                ntree = best_ntree)


confusion_table<-cbind(rownames(best_fit$finalModel$confusion), best_fit$finalModel$confusion)
write.table(confusion_table,  paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_train_confusion_matrix.txt", sep=''), row.names=F, sep="\t", na="", quote=F)

save(best_fit,file = paste("covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_model.RData", sep=''))

library(pROC)
library(ROCR)
pred<-best_fit$pred
pp<-c()
for(i in 1:3){
  if(i==1){
    aa<-pred[grep(paste("Rep", i, sep=''), pred$Resample), c("rowIndex", "Y", "obs")]
    pp<-aa
  }else{
    aa<-pred[grep(paste("Rep", i, sep=''), pred$Resample), c("rowIndex", "Y")]
    pp<-merge(aa,pp, by="rowIndex")
  }
}

pp$prob<-apply(pp[,2:4],1,mean)
rf_t_prob <- pp[order(pp$rowIndex, decreasing = F), c("rowIndex", "obs", "prob")]
rownames(rf_t_prob)<-rownames(genes_mat)
roc_curve<-roc(rf_t_prob$obs, rf_t_prob$prob)

pdf( paste( "covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_train_ROC.pdf", sep=""), width=6, height=6)
plot(roc_curve, type="l", xlab="False Positive Rate",ylab="Ture Positive Rate", col="black")
legend(0.5 ,0.2, paste("AUC = ",round(as.numeric(gsub("Area under the curve: ", "", roc_curve$auc)),3), "", sep=''),
       col=c("black", "blue"), cex=0.8, pch=1, lty=1)
dev.off()




#apply to external data set
new_data<-read.delim(paste("covid_proteomics_adjusted_data.txt", sep=''), check.names = F)
colnames(new_data)<-gsub("_Run.*", "", colnames(new_data))
rownames(new_data)<-new_data$Accession

sample_test<-read.delim("covid_proteomics_samples_sum.txt")
sample_test<-sample_test[sample_test$Set=="test"&!is.na(sample_test$Composite_Outcome),]
sample_test$Study_ID<-sample_test$Sample
new_data<-new_data[(new_data$Accession%in%colnames(genes_mat)), as.character(sample_test$Study_ID)]
new_data2<-t(new_data)
new_data2<-new_data2[, colnames(genes_mat)[-ncol(genes_mat)]]


rf_v_prob <- predict(best_fit, newdata=new_data2, type="prob")
rf_v_prob$Sample=rownames(rf_v_prob)
rf_v_prob<-merge(sample_test, rf_v_prob[, c("Sample", "Y")])
colnames(rf_v_prob)[ncol(rf_v_prob)]="prob"

rownames(rf_v_prob)<-rf_v_prob$Sample
roc_curve<-roc(rf_v_prob$Composite_Outcome, rf_v_prob$prob)

pdf( paste( "covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_test_ROC.pdf", sep=""), width=6, height=6)
plot(roc_curve, type="l", xlab="False Positive Rate",ylab="Ture Positive Rate", col="black")
legend(0.5 ,0.2, paste("AUC = ",round(as.numeric(gsub("Area under the curve: ", "", roc_curve$auc)),3), "", sep=''),
       col=c("black", "blue"), cex=0.8, pch=1, lty=1)
dev.off()



sample_info<-read.delim("covid_proteomics_samples_sum.txt")
sample_info<-sample_info[sample_info$Set=="train"&!is.na(sample_info$Composite_Outcome),]
sample_info$Study_ID<-sample_info$Sample

rf_t_prob2<-merge(data.frame(Study_ID=rownames(rf_t_prob), rf_t_prob), sample_info, by="Study_ID")
rf_t_prob2<-rf_t_prob2[order(rf_t_prob2$prob, decreasing = T),]

rf_v_prob<-rf_v_prob[order(rf_v_prob$prob, decreasing = T),]

library(xlsx)
write.xlsx(rf_t_prob2, file=paste( "covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_prediction.xlsx", sep=""),
           sheetName="training", append=F)
write.xlsx(rf_v_prob, file=paste( "covid_proteomics_Composite_Outcome_Y_vs_N_unadj_limma_0.05p_feature_sel_RF_prediction.xlsx", sep=""),
           sheetName="testing", append=TRUE)