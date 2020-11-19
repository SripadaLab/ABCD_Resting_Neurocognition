setwd("/net/pepper/ABCD/CIFTI/Scripts/rest_neurocognition/")
source("setup.R")

library(lme4)

# read phenotypic data file
pheno = read.csv(file.path(DataDir,"ABCD_rest.csv"),na.strings=c("NA","NaN"," "))

# read PCA expressions
expressions = read.csv(file.path(ResultsDir,"ABCD_rest_CIFTI_expressions.csv"),na.strings=c("NA","NaN"," "))

# merge
dat = merge(pheno,expressions,by="subjectkey")

controls = c('Gender', 'Age','I(Age^2)','C(RaceEthnicity)','fd', 'I(fd^2)')

# and multilevel structure
nested_str = '(1| site_id_l/rel_family_id)'

n_components=250
brain_formula='A001'
for (i in 2:n_components) {
     icomp <- sprintf("%03d",i) # fix to 3 characters 
     brain_formula = paste(brain_formula, ' + ', 'A', icomp, sep = "")
}

outcomes = c("G_lavaan","S1_lavaan","S2_lavaan","S3_lavaan",
  "nihtbx_cardsort_uncorrected","nihtbx_flanker_uncorrected",
  "nihtbx_list_uncorrected","nihtbx_pattern_uncorrected",
  "nihtbx_picture_uncorrected","nihtbx_picvocab_uncorrected",
  "nihtbx_reading_uncorrected","pea_ravlt_sd_tc",
  "pea_ravlt_ld_tc","pea_wiscv_trs","lmt_scr_num_correct")

betas = matrix(NA,nrow=250,ncol=length(outcomes))
ts = matrix(NA,nrow=250,ncol=length(outcomes))

for (i in 1:length(outcomes)) {
  dependent = outcomes[i]
  fmla = paste(dependent, ' ~ ', brain_formula, ' + ', 
    paste(controls, collapse=' + '), ' + ', 
    nested_str, sep='')
  
  model = lmer(formula=fmla, data=dat, na.action='na.exclude')
  
  betas[,i] = coef(summary(model))[2:251,1]
  ts[,i] = coef(summary(model))[2:251,3]
}

write.table(betas,file.path(ResultsDir,"mm_betas.csv"),row.names=F,col.names=F,quote=F,na="NaN")

dts = as.data.frame(ts)
names(dts) = outcomes
dts$Component = 1:250
dts = dts[,c("Component",outcomes)]
write.csv(dts,file.path(ResultsDir,"ABCD_component_ts.csv"),row.names=F,quote=F)
