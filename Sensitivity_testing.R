################################
library(dplyr)

##Using olink MISC and healthy to compare 

olink_data=read.table("/Data/COVID_NPX.csv", sep=",", header=TRUE)

head(olink_data)
rownames(olink_data)=olink_data$SampleID
olink_data$SampleID=NULL
#Separating timepoint 1 from timepoint 2
tp1=olink_data[which(olink_data$Timepoint==1),]

head(tp1)
tp2=olink_data[which(olink_data$Timepoint==2),]
head(tp2)

tp1_pca=select(tp1, -(1:4))
head(tp1_pca)

#Looking at increasing healthy values and what the difference in DE is 

tp1_1.1=tp1_pca
tp1_1.1[which(tp1$DiseaseCat == 4),]=1.1*tp1_1.1[which(tp1$DiseaseCat==4),]

tp1_1.2=tp1_pca
tp1_1.2[which(tp1$DiseaseCat == 4),]=1.2*tp1_1.2[which(tp1$DiseaseCat==4),]

tp1_1.3=tp1_pca
tp1_1.3[which(tp1$DiseaseCat == 4),]=1.3*tp1_1.3[which(tp1$DiseaseCat==4),]

tp1_1.4=tp1_pca
tp1_1.4[which(tp1$DiseaseCat == 4),]=1.4*tp1_1.4[which(tp1$DiseaseCat==4),]

de_groups2 = function(data, group1,group2){
  diff_exp_matrix_fc=apply(data,MARGIN=2,function(X)
    mean(X[which(tp1$DiseaseCat==group1)])-mean(X[which(tp1$DiseaseCat==group2)]))
  diff_exp_matrix_pval=apply(data,MARGIN=2,function(X)
    p.adjust(t.test(X[which(tp1$DiseaseCat==group1)],X[which(tp1$DiseaseCat==group2)])$p.value,method="BH",n=length(5:(ncol(tp1)-1))))
  #filtering based on pval 
  Pval_Threshold1 = 0.05
  FC_Threshold = 2
  Sig_text = names(diff_exp_matrix_pval)[which(diff_exp_matrix_pval<= Pval_Threshold1)]
  #Filtering based on FC 
  FC_text=names(diff_exp_matrix_fc)[which(abs(diff_exp_matrix_fc)>= FC_Threshold)]
  DE_Genes=intersect(Sig_text, FC_text)
  de_fc=diff_exp_matrix_fc[DE_Genes]
  de_pval=diff_exp_matrix_pval[DE_Genes]
  final_mat=data.frame(
    "prot.name"= DE_Genes,
    "fc"= de_fc,
    "Pval"= de_pval
  )
  return(final_mat)
  
}

tp1$DiseaseCategory=tp1$DiseaseCat
tp1$DiseaseCategory[which(tp1$DiseaseCategory == "1")] = "MISC"
tp1$DiseaseCategory[which(tp1$DiseaseCategory == "2")] = "Minimal"
tp1$DiseaseCategory[which(tp1$DiseaseCategory == "3")] = "Severe"
tp1$DiseaseCategory[which(tp1$DiseaseCategory == "4")] = "Healthy"
head(tp1$DiseaseCategory)

##MISC vs Healthy 

DE_1.4=de_groups2(tp1_pca, "1", "4")
DE_1.4_1.1=de_groups2(tp1_1.1, "1", "4")
DE_1.4_1.2=de_groups2(tp1_1.2, "1", "4")
DE_1.4_1.3=de_groups2(tp1_1.3, "1", "4")
DE_1.4_1.4=de_groups2(tp1_1.4, "1", "4")

#Minimal vs Healthy 
DE_2.4=de_groups2(tp1_pca, "2", "4")
DE_2.4_1.1=de_groups2(tp1_1.1, "2", "4")
DE_2.4_1.2=de_groups2(tp1_1.2, "2", "4")
DE_2.4_1.3=de_groups2(tp1_1.3, "2", "4")
DE_2.4_1.4=de_groups2(tp1_1.4, "2", "4")


#Severe vs Healthy 
DE_3.4=de_groups2(tp1_pca, "3", "4")
DE_3.4_1.1=de_groups2(tp1_1.1, "3", "4")
DE_3.4_1.2=de_groups2(tp1_1.2, "3", "4")
DE_3.4_1.3=de_groups2(tp1_1.3, "3", "4")
DE_3.4_1.4=de_groups2(tp1_1.4, "3", "4")


testing_de=function(de1, de2, slopegroup){
  all_genes=unique(c(de1$prot.name, de2$prot.name, slopegroup))
  out=data.frame(genes=all_genes)
  out$group=rep("", nrow(out))
  both=all_genes[which(de1$prot.name %in% de2$prot.name)]
  for(i in 1:length(both)){
    if(de1$fc[which(de1$prot.name == both[i])] > 0 & de2$fc[which(de2$prot.name == both[i])] > 0 |
       de1$fc[which(de1$prot.name == both[i])] < 0 & de2$fc[which(de2$prot.name == both[i])] < 0 ){
      out$group[which(out$genes == both[i])] = "TP"
    }else{out$group[which(out$genes == both[i])] = 'FP'}
  }
  other_genes=all_genes[-which(all_genes %in% both)]
  for( o in 1:length(other_genes)){
    if(other_genes[o] %in% de2$prot.name){out$group[which(out$genes == other_genes[o])]='FP'}
    else{ out$group[which(out$genes == other_genes[o])]='FN'}
  }
  out$slope=slopegroup
  return(out)
}


diff_de1.4_1.1=testing_de(DE_1.4, DE_1.4_1.1, 'slope_1.1')
diff_de1.4_1.2=testing_de(DE_1.4, DE_1.4_1.2, 'slope_1.2')
diff_de1.4_1.3=testing_de(DE_1.4, DE_1.4_1.3, 'slope_1.3')
diff_de1.4_1.4=testing_de(DE_1.4, DE_1.4_1.4, 'slope_1.4')

de_analysis=rbind(diff_de1.4_1.1, diff_de1.4_1.2, diff_de1.4_1.3, diff_de1.4_1.4)

length(which(diff_de1.4_1.4$group %in% c("FP")))/197
length(which(diff_de1.4_1.4$group %in% c('FN')))/212
length(which(diff_de1.4_1.4$group %in% c('TP')))/212

length(which(diff_de1.4_1.3$group %in% c("FP")))/201
length(which(diff_de1.4_1.3$group %in% c("FN")))/212

length(which(diff_de1.4_1.2$group %in% c("FP", 'FN')))
length(which(diff_de1.4_1.1$group %in% c("FP", 'FN')))

31/247

library(ggplot2)
ggplot(de_analysis, aes(x=slope)) + geom_bar(aes(fill=group)) + theme_classic() + ggtitle('DEA MISC vs. Healthy - Positive Slope') + ylab('Protein Count') + scale_fill_brewer(palette = 'YlOrRd',direction = -1 ) + 
  theme(text=element_text(size=25)) + xlab('')

diff_de2.4_1.1=testing_de(DE_2.4, DE_2.4_1.1, 'slope_1.1')
diff_de2.4_1.2=testing_de(DE_2.4, DE_2.4_1.2, 'slope_1.2')
diff_de2.4_1.3=testing_de(DE_2.4, DE_2.4_1.3, 'slope_1.3')
diff_de2.4_1.4=testing_de(DE_2.4, DE_2.4_1.4, 'slope_1.4')

length(which(diff_de2.4_1.1$group %in% c("FP", 'FN')))
length(which(diff_de2.4_1.2$group %in% c("FP", 'FN')))
length(which(diff_de2.4_1.3$group %in% c("FP")))/60
length(which(diff_de2.4_1.3$group %in% c("FN")))/75
length(which(diff_de2.4_1.3$group %in% c("TP")))/75

length(which(diff_de2.4_1.4$group %in% c("FP")))/64
length(which(diff_de2.4_1.4$group %in% c('FN')))/75
length(which(diff_de2.4_1.4$group %in% c('TP')))/75

62/101

de_analysis2=rbind(diff_de2.4_1.1, diff_de2.4_1.2, diff_de2.4_1.3, diff_de2.4_1.4)



ggplot(de_analysis2, aes(x=slope)) + geom_bar(aes(fill=group)) + theme_classic() + ggtitle('DEA Minimal COVID vs. Healthy - Positive Slope') + ylab('Protein Count') + 
  scale_fill_brewer(palette = 'YlOrRd',direction = -1 ) + 
  theme(text=element_text(size=25), legend.position = 'right') + xlab('') 


#################################################################################
#Looking at negative slope 

tp1_1.1=tp1_pca
tp1_1.1[which(tp1$DiseaseCat == 4),]=tp1_1.1[which(tp1$DiseaseCat==4),]*0.9

tp1_1.2=tp1_pca
tp1_1.2[which(tp1$DiseaseCat == 4),]=tp1_1.2[which(tp1$DiseaseCat==4),]*0.8

tp1_1.3=tp1_pca
tp1_1.3[which(tp1$DiseaseCat == 4),]=tp1_1.3[which(tp1$DiseaseCat==4),]*0.7

tp1_1.4=tp1_pca
tp1_1.4[which(tp1$DiseaseCat == 4),]=tp1_1.4[which(tp1$DiseaseCat==4),]*0.6

DE_1.4=de_groups2(tp1_pca, "1", "4")
DE_1.4_0.9=de_groups2(tp1_1.1, "1", "4")
DE_1.4_0.8=de_groups2(tp1_1.2, "1", "4")
DE_1.4_0.7=de_groups2(tp1_1.3, "1", "4")
DE_1.4_0.6=de_groups2(tp1_1.4, "1", "4")

DE_2.4=de_groups2(tp1_pca, "2", "4")
DE_2.4_0.9=de_groups2(tp1_1.1, "2", "4")
DE_2.4_0.8=de_groups2(tp1_1.2, "2", "4")
DE_2.4_0.7=de_groups2(tp1_1.3, "2", "4")
DE_2.4_0.6=de_groups2(tp1_1.4, "2", "4")


DE_3.4=de_groups2(tp1_pca, "3", "4")
DE_3.4_1.1=de_groups2(tp1_1.1, "3", "4")
DE_3.4_1.2=de_groups2(tp1_1.2, "3", "4")
DE_3.4_1.3=de_groups2(tp1_1.3, "3", "4")
DE_3.4_1.4=de_groups2(tp1_1.4, "3", "4")


diff_de1.4_0.9=testing_de(DE_1.4, DE_1.4_0.9, 'slope_0.9')
diff_de1.4_0.8=testing_de(DE_1.4, DE_1.4_0.8, 'slope_0.8')
diff_de1.4_0.7=testing_de(DE_1.4, DE_1.4_0.7, 'slope_0.7')
diff_de1.4_0.6=testing_de(DE_1.4, DE_1.4_0.6, 'slope_0.6')

length(which(diff_de1.4_0.9$group %in% c("FP", 'FN')))
length(which(diff_de1.4_0.8$group %in% c("FP", 'FN')))

length(which(diff_de1.4_0.7$group %in% c("FP")))/227
length(which(diff_de1.4_0.7$group %in% c("FN")))/212
length(which(diff_de1.4_0.7$group %in% c("TP")))/212

length(which(diff_de1.4_0.6$group %in% c("FP")))/234
length(which(diff_de1.4_0.6$group %in% c('FN')))/212
length(which(diff_de1.4_0.6$group %in% c('TP')))/212

de_analysis4=rbind(diff_de1.4_0.9, diff_de1.4_0.8, diff_de1.4_0.7, diff_de1.4_0.6)

dark_palette <- brewer.pal(9, "YlGnBu")[3:9]  # Adjust the indices to select the darker colors


ggplot(de_analysis4, aes(x=slope)) + 
  geom_bar(aes(fill=group)) + 
  theme_classic() + 
  ggtitle('DEA MISC vs. Healthy - Negative Slope') + 
  ylab('Protein Count') + 
  scale_fill_manual(values = dark_palette) + 
  theme(text = element_text(size=25), legend.position = 'right') + 
  xlab('')

diff_de2.4_0.9=testing_de(DE_2.4, DE_2.4_0.9, 'slope_0.9')
diff_de2.4_0.8=testing_de(DE_2.4, DE_2.4_0.8, 'slope_0.8')
diff_de2.4_0.7=testing_de(DE_2.4, DE_2.4_0.7, 'slope_0.7')
diff_de2.4_0.6=testing_de(DE_2.4, DE_2.4_0.6, 'slope_0.6')

length(which(diff_de2.4_0.9$group %in% c("FP", 'FN')))
length(which(diff_de2.4_0.8$group %in% c("FP", 'FN')))
length(which(diff_de2.4_0.7$group %in% c("FP", 'FN')))
length(which(diff_de2.4_0.6$group %in% c("FP")))



de_analysis5=rbind(diff_de2.4_0.9, diff_de2.4_0.8, diff_de2.4_0.7, diff_de2.4_0.6)

ggplot(de_analysis5, aes(x=slope)) + geom_bar(aes(fill=group)) + theme_classic() + ggtitle('DEA Minimal COVID vs. Healthy - Negative Slope') + ylab('Protein Count') + 
  scale_fill_manual(values = dark_palette) + 
  theme(text = element_text(size=25), legend.position = 'bottom') + 
  xlab('')
