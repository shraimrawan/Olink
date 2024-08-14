#########################################
#Validation cohort 

npx_val=read.delim("/Data/Validation_NPX.csv", sep=",", header=T)

models_allprot=read.delim("/Data/finalmodelsprotein.csv", sep=",", row.names = "X")

meta=read.delim("/Data/ValidationCohort_Decoder.csv", sep=",", header=T)


models_allprot$slope.std=round(models_allprot$slope.std, 2)

lowcimodels_orig=models_allprot[which(abs(models_allprot$slope.std) <= 0.3),]


colnames(npx_val)[-which(colnames(npx_val) %in% rownames(models_allprot))]

#Making Olink protein names consistent across datasets

colnames(npx_val)[which(colnames(npx_val) == "MCP.3")]="CCL7"
colnames(npx_val)[which(colnames(npx_val) == "MCP.1")]="CCL2"
colnames(npx_val)[which(colnames(npx_val) == "MCP.4")]="CCL13" 
colnames(npx_val)[which(colnames(npx_val) == "TRAIL")]="TNFSF10" 
colnames(npx_val)[which(colnames(npx_val) == "CAIX")]="CA9" 
colnames(npx_val)[which(colnames(npx_val) == "Gal.9")]="LGALS9"
colnames(npx_val)[which(colnames(npx_val) == "PDGF.subunit.B")]="PDGFRB"
colnames(npx_val)[which(colnames(npx_val) == "MCP.2")]="CCL8"
colnames(npx_val)[which(colnames(npx_val) == "PD.L1")]="CD274"
colnames(npx_val)[which(colnames(npx_val) == "HO.1")]="HMOX1"
colnames(npx_val)[which(colnames(npx_val) == "MIC.A.B")]="MICB_MICA"
colnames(npx_val)[which(colnames(npx_val) == "IFN.gamma")]="IFNG"
colnames(npx_val)[which(colnames(npx_val) == "CSF.1")]="CSF1"
colnames(npx_val)[which(colnames(npx_val) == "TIE2")]="TEK"
colnames(npx_val)[which(colnames(npx_val) == "IL12")]="IL12A_IL12B"

valid_prot=colnames(npx_val)[which(colnames(npx_val) %in% rownames(models_allprot))]


plasma_meta=meta[which(meta$Substrate == "Plasma"),]

plasma_meta2=plasma_meta[-which(plasma_meta$SampleID == "Sample control"),]

serum_meta=meta[which(meta$Substrate == "Serum"),]

serum_meta=serum_meta[order(serum_meta$SampleID),]

plasma_meta2=plasma_meta2[order(plasma_meta2$SampleID),]

serum_npx=npx_val[which(npx_val$SampleID %in% serum_meta$OlinkID),]
plasma_npx=npx_val[which(npx_val$SampleID %in% plasma_meta2$OlinkID),]

modeled_prot=models_allprot[which(rownames(models_allprot) %in% colnames(npx_val)),]

rownames(plasma_npx)=plasma_npx$SampleID

plasma_npx2=plasma_npx[order(match(plasma_npx$SampleID, plasma_meta2$OlinkID)),]
serum_npx2=serum_npx[order(match(serum_npx$SampleID, serum_meta$OlinkID)),]

proteins_modeled=rownames(models_allprot)
proteins_allpanel=colnames(serum_npx2)[2:93]

corr_pearson=data.frame()

for (i in 1:length(proteins_allpanel)){
  cor_cof=cor(plasma_npx2[,proteins_allpanel[i]], serum_npx2[,proteins_allpanel[i]], method="pearson", use="complete.obs")
  corr_pearson=rbind(corr_pearson, cor_cof)
}

row.names(corr_pearson)=proteins_allpanel
colnames(corr_pearson)="correlation_pearson"

length(which(abs(corr_pearson$correlation_pearson) >= 0.5))

valid_linear=linear_model(proteins_allpanel)
valid_linear$R.squared=round(valid_linear$R.squared, 1)
valid_linear$tier='Tier1'

length(which(valid_linear$R.squared >= 0.5))

valid_linear_rm2=linear_rmoutlier(2, proteins_allpanel)
valid_linear_rm2$R.squared=round(valid_linear_rm2$R.squared,1 )

length(which(valid_linear_rm2$R.squared >= 0.5))
valid_linear_rm2$tier='Tier2'


valid_linear_rm3=linear_rmoutlier(3, proteins_allpanel)
valid_linear_rm3$R.squared=round(valid_linear_rm3$R.squared,1)

length(which(valid_linear_rm3$R.squared >= 0.5))
valid_linear_rm3$tier='Tier3'


a=rownames(valid_linear)[which(valid_linear$R.squared >= 0.5)]
b=rownames(valid_linear_rm2)[which(valid_linear_rm2$R.squared >= 0.5)]
c=rownames(valid_linear_rm3)[which(valid_linear_rm3$R.squared >= 0.5)]

comb=c(a,b,c)
com_all=(unique(comb))


b_prot=b[-which(b %in% a)]
c_prot=c[-which(c %in% c(b,a))]

valid_models=valid_linear[which(rownames(valid_linear) %in% a),]
valid_models=rbind(valid_models, valid_linear_rm2[which(rownames(valid_linear_rm2) %in% b_prot),])
valid_models=rbind(valid_models, valid_linear_rm3[which(rownames(valid_linear_rm3) %in% c_prot),])


rownames(valid_models)[which(rownames(valid_models) == 'MUC.16')]='MUC16'
rownames(valid_models)[which(rownames(valid_models) == 'CASP.8')]='CASP8'
rownames(valid_models)[which(rownames(valid_models) == 'Gal.1')]='LGALS1'
rownames(valid_models)[which(rownames(valid_models) == 'CD40.L')]='CD40LG'
rownames(valid_models)[which(rownames(valid_models) == 'PD.L2')]='PDCD1LG2'
rownames(valid_models)[which(rownames(valid_models) == 'TWEAK')]='TNFSF12'

lowcivalid=valid_models[-which(rownames(valid_models) %in% rownames(lowcimodels)),]


##########
#Looking at proteins with highCI and non-modelable 

modeled_prots=c(rownames(lowcimodels), rownames(lowcivalid))

colnames(lowcivalid)[5]='lm_tier'

final_factors=rbind(lowcimodels, lowcivalid)

nomodel=slopes_all[-which(rownames(slopes_all) %in% modeled_prots),]

a=as.data.frame(table(modeled_prots))

#Supplement table from paper 
paper_prots=read.csv("~/Library/CloudStorage/Box-Box/Teachey_Lab/Projects/Olinking/Data/Olink_serumplasma_paper.csv", sep = ",")

nomodel_paper=rownames(nomodel)[which(rownames(nomodel) %in% paper_prots$proteins)]

model_paper=modeled_prots[which(modeled_prots %in% paper_prots$proteins)]

model_paper_factors=final_factors[model_paper,]

corr_pearson_nomodelpaper=corr_pearson[nomodel_paper, ,drop=F]

corr_pearson_nomodelpaper$correlation_pearson=round(corr_pearson_nomodelpaper$correlation_pearson, 2)

corr_pearson_nomodelpaper_high=corr_pearson_nomodelpaper[corr_pearson_nomodelpaper$correlation_pearson >= 0.5,, drop=F]

corr_pearson_nomodelpaper_high=corr_pearson_nomodelpaper[corr_pearson_nomodelpaper$correlation_pearson >= 0.8,, drop=F]


t_test_results <- data.frame(Protein = nomodel_paper)

t_test_results$p_value <- sapply(nomodel_paper, function(p) {
  t.test(serum_dirt[,p], plasma_dirt[,p], paired = TRUE)$p.value
})

which(t_test_results$p_value <= 0.05)


t_test_results2 <- data.frame(Protein = all_prot)

t_test_results2$p_value <- sapply(all_prot, function(p) {
  t.test(serum_dirt[,p], plasma_dirt[,p], paired = TRUE)$p.value
})

which(t_test_results2$p_value <= 0.05)

models_allprot$cigroup=rep('<= 0.3', nrow(models_allprot))
models_allprot$cigroup[which(abs(models_allprot$slope.std) >= 0.3)]='> 0.3'

ggplot(models_allprot, aes(x=lm_tier)) + geom_bar(aes(fill=cigroup)) + theme_classic()  + ylab('Protein Count') + 
  theme(text = element_text(size=25)) + scale_fill_manual(values = c('darkgreen', 'grey')) + 
  guides(fill=guide_legend(title="Slope CI")) + xlab('')

