####################
##Reading in data 
rm(list = ls())

serum_dirt=read.delim("/Data/normalized_rosetta_serum_NPX.csv", sep=",", header=T)

plasma_dirt=read.delim("/Data/normalized_rosetta_plasma_NPX.csv", sep=",", header=T)

decoder=read.delim("/Data/InitialCohort_Decoder.csv", sep=",", header=T)

decoder=decoder[1:19,]

all_prot=colnames(plasma_dirt)[2:1464]

####################
##Functions to plot data 
library(ggpmisc)

ploting_olink=function(protein){
  df=data.frame(x=plasma_dirt[,protein], y=serum_dirt[,protein])
  ggplot(df, aes(x=x, y=y)) + geom_point() + theme_bw() + ylab(paste0("plasma_", protein)) + xlab(paste0("serum_", protein)) + geom_smooth(method=lm) + stat_poly_eq() + ggtitle('DIRT Cohort')
}

####################
##Calculating pearson correlation 

corr_pearson=data.frame()

for (i in 1:length(all_prot)){
  cor_cof=cor(plasma_dirt[,all_prot[i]], serum_dirt[,all_prot[i]], method="pearson", use="complete.obs")
  corr_pearson=rbind(corr_pearson, cor_cof)
}

row.names(corr_pearson)=all_prot
colnames(corr_pearson)="correlation_pearson"

length(which(corr_pearson$correlation_pearson <= -0.5))

####################
##Linear modeling 

linear_model=function(proteins){
  final_df=data.frame()
  for(i in 1: length(proteins)){
    df_mod=data.frame(x=serum_dirt[,proteins[i]], y=plasma_dirt[, proteins[i]])
    linear_eq=lm(x ~ y, data = df_mod)
    out=summary(linear_eq)
    std=out$coefficients[2,2]*2
    b=linear_eq$coefficients[1]
    slope=linear_eq$coefficients[2]
    rsquared=summary(linear_eq)$r.squared
    add=c(b, slope, rsquared, std)
    final_df=rbind(final_df, add)
    colnames(final_df)=c("Intercept", "Slope", "R.squared", "slope.std")
  }
  rownames(final_df)=proteins
  
  return(final_df)
}


##Linear modeling removing outliers 

linear_rmoutlier=function(outliers_rm, proteins){
  final_df=data.frame()
  for(i in 1: length(proteins)){
    df_mod=data.frame(x=serum_dirt[,proteins[i]], y=plasma_dirt[, proteins[i]])
    linear_eq=lm(x ~ y, data = df_mod)
    cooksd=cooks.distance(linear_eq)
    samplesize=nrow(df_mod)
    influential=as.numeric(names(sort(cooksd, decreasing = TRUE)[1:outliers_rm]))
    df2=df_mod[-influential,]
    linear_eq2=lm(x ~ y, data=df2)
    out=summary(linear_eq2)
    std=out$coefficients[2,2]*2
    b=linear_eq2$coefficients[1]
    slope=linear_eq2$coefficients[2]
    rsquared=summary(linear_eq2)$r.squared
    add=c(b, slope, rsquared, std)
    final_df=rbind(final_df, add)
    colnames(final_df)=c("Intercept", "Slope", "R.squared", "slope.std")
  }
  rownames(final_df)=proteins
  return(final_df)
}

slopes_all=linear_model(all_prot)

slopes_all$R.squared=round(slopes_all$R.squared, 1)
length(which(slopes_all$R.squared >= 0.5))

outliers2=linear_rmoutlier(2, all_prot)
outliers2$R.squared=round(outliers2$R.squared, 1)


outliers3=linear_rmoutlier(3, all_prot)

outliers3$R.squared=round(outliers3$R.squared, 1)

tier1=rownames(slopes_all)[which(abs(slopes_all$R.squared) >= 0.5)]

tier2=rownames(outliers2)[which(abs(outliers2$R.squared) >= 0.5)]

tier2=tier2[-which(tier2 %in% tier1)]

tier3=rownames(outliers3)[which(abs(outliers3$R.squared) >= 0.5)]

tier3=tier3[-which(tier3 %in% c(tier1, tier2))]

linear_tier1=slopes_all[which(slopes_all$R.squared >= 0.5),]

low_liner=slopes_all[which(slopes_all$R.squared < 0.5),]

linear_out2r=outliers2[which(outliers2$R.squared >= 0.5 & rownames(outliers2) %in% rownames(low_liner)),]

linear_out3r=outliers3[which(outliers3$R.squared >= 0.5 & rownames(outliers3) %in% rownames(low_liner)),]

a=which(rownames(linear_out3r) %in% rownames(linear_out2r))

linear_out3r=linear_out3r[-a,]

##########

#Looking at proteins that were only modeled in tier 2 or tier 3 
slopes_all=linear_model(all_prot)
outliers2=linear_rmoutlier(2, all_prot)
outliers3=linear_rmoutlier(3, all_prot)

lower_models=slopes_all[c(tier2, tier3),c('R.squared', 'slope.std')]

lower_models2=merge(lower_models, outliers2[,3:4], by=0, all.x = T)

lower_models3=merge(lower_models2, outliers3[,3:4], by.x = 'Row.names', by.y=0, all.x = T)

colnames(lower_models3)=c('Proteins', 'Rsquared.Tier1', 'STD.Tier1', 'Rsquared.Tier2', 'STD.Tier2', 'Rsquared.Tier3', 'STD.Tier3')

tier_2=lower_models3[which(lower_models3$rs.tier2 >= 0.5),]

tier_3=lower_models3[which(lower_models3$Proteins %in% tier3),]

library(reshape2)

summary=melt(tier_3)

summary.rs=summary[which(grepl('Rsquared', summary$variable)),]
summary.std=summary[which(grepl('STD', summary$variable)),]
summary.rs$Proteins=as.factor(summary.rs$Proteins)
summary.std$Proteins=as.factor(summary.std$Proteins)

library(RColorBrewer)
library(ggplot2)

num_colors <- length(unique(summary.rs$Proteins))

# Create a custom color palette
custom_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)

ggplot(summary.rs, aes(x=variable, y=value, group=Proteins, color=Proteins)) + geom_point() + geom_line() +theme_classic() + ylab('R-squared') + xlab(' Linear Model Tiers') + scale_color_manual(values = custom_palette) + 
  theme(text = element_text(size=15))

ggplot(summary.std, aes(x=variable, y=value, group=Proteins, color=Proteins)) + geom_point() + geom_line() +theme_classic() + ylab('R-squared (rs)') + xlab(' Linear Model Tiers')


##################################
## Slope and CI of modelable proteins 

linear_tier1$lm_tier='Tier1'
linear_out2r$lm_tier='Tier2'
linear_out3r$lm_tier='Tier3'

library(ggplot2)

models_allprot=rbind(linear_tier1, linear_out2r, linear_out3r)

lowcimodels=models_allprot[which(models_allprot$slope.std >= -0.3 & models_allprot$slope.std <= 0.3), ]

highcimodels=models_allprot[which(models_allprot$slope.std < -0.3 | models_allprot$slope.std > 0.3), ]

length(which(highcimodels$lm_tier == 'Tier1' ))
length(which(highcimodels$lm_tier == 'Tier2' ))
length(which(highcimodels$lm_tier == 'Tier3' ))

ggplot(highcimodels, aes(x=lm_tier)) + geom_bar(stat = 'count', aes(fill=lm_tier)) + theme_classic()


###Slope CI graph
purple_palette <- brewer.pal(8, "Purples")

# Create the histogram with custom fill colors
ggplot(models_allprot, aes(x=slope.std)) + 
  geom_histogram(binwidth = 0.1, fill='seagreen3') + # Use a specific shade of purple
  ylab('Count of Proteins') + 
  xlab('Slope CI') + 
  theme_classic() + 
  theme(text = element_text(size=15))

####Slope value graph
ggplot(models_allprot, aes(x=Slope)) + 
  geom_histogram(binwidth = 0.2) + ylab('Count of Proteins') + xlab('Slope') + theme_classic() + theme(text = element_text(size=12))

ggplot(models_allprot, aes(x=Slope)) + 
  geom_histogram(binwidth = 0.1, fill=purple_palette[5]) + # Use a specific shade of purple
  ylab('Count of Proteins') + 
  xlab('Slope') + 
  theme_classic() + 
  theme(text = element_text(size=15))


length(which(models_allprot$slope.std > 0 & models_allprot$slope.std <= 0.3))

length(which(models_allprot$slope.std < 0 ))


models_allprot$Slope.CI=rep('LowCI', nrow(models_allprot))

models_allprot$Slope.CI[which(models_allprot$slope.std < -0.3 | models_allprot$slope.std > 0.3)]='HighCI'

table(models_allprot$Slope.CI)

ggplot(models_allprot, aes(x=lm_tier)) + geom_bar(position = 'fill', aes(fill=Slope.CI)) + ylab('Proteins (%)') + theme_classic() + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = c('peachpuff', 'lightgreen')) + 
  theme(text=element_text(size=15)) + xlab('')


