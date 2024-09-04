source('abcd_functions.R')
library(lavaan)

ABCDDataDir = '../TextFiles/'

# 1. "nihtbx_picvocab_uncorrected", 
# 2. "nihtbx_flanker_uncorrected",
# 3. "nihtbx_list_uncorrected", 
# 4. "nihtbx_cardsort_uncorrected",
# 5. "nihtbx_pattern_uncorrected",
# 6. "nihtbx_picture_uncorrected",
# 7. "nihtbx_reading_uncorrected", 
# 8. "lmt_scr_num_correct",
# 9.  "pea_ravlt_sd_tc", 
# 10. "pea_ravlt_ld_tc"
# 11. "pea_wiscv_trs")

nih = read.abcd(file.path(ABCDDataDir,"abcd_tbss01.txt"))
nih = nih[,c("subjectkey","eventname","nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected",
             "nihtbx_list_uncorrected","nihtbx_cardsort_uncorrected","nihtbx_pattern_uncorrected",
             "nihtbx_picture_uncorrected","nihtbx_reading_uncorrected",
             "nihtbx_cryst_uncorrected","nihtbx_fluidcomp_uncorrected")]

ps = read.abcd(file.path(ABCDDataDir,"abcd_ps01.txt"))
ps$pea_ravlt_sd_tc = rowSums(ps[,names(ps)[grepl("pea_ravlt_sd.*_tc",names(ps))]])
ps$pea_ravlt_ld_tc = ps$pea_ravlt_ld_trial_vii_tc
ps = ps[,c("subjectkey","eventname","pea_ravlt_sd_tc","pea_ravlt_ld_tc","pea_wiscv_trs")]

lmt = read.abcd(file.path(ABCDDataDir,"lmtp201.txt"))
lmt = lmt[,c("subjectkey","eventname","lmt_scr_num_correct")]

site = read.abcd(file.path(ABCDDataDir,"abcd_lt01.txt"))
site = site[,c("subjectkey","eventname","site_id_l")]

data = multi.merge(site,nih,ps,lmt,by=c("subjectkey","eventname"))
data = data[data$eventname=="baseline_year_1_arm_1",]

bad = apply(is.na(data[,c(4,5,6,7,8,9,10,13,14,15,16)]),1,any)

data = data[!bad,]

u = sort(unique(data$site_id_l))
nFold = length(u)
cols = c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_list_uncorrected",
         "nihtbx_cardsort_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected",
         "nihtbx_reading_uncorrected","pea_ravlt_sd_tc","pea_ravlt_ld_tc","pea_wiscv_trs",
         "lmt_scr_num_correct")

bifactor = '
G =~ nihtbx_picvocab_uncorrected + nihtbx_flanker_uncorrected + nihtbx_list_uncorrected + nihtbx_cardsort_uncorrected + 
    nihtbx_pattern_uncorrected + nihtbx_picture_uncorrected + nihtbx_reading_uncorrected + lmt_scr_num_correct  +  pea_ravlt_sd_tc  +  pea_ravlt_ld_tc +  pea_wiscv_trs
S1 =~ nihtbx_picvocab_uncorrected + nihtbx_list_uncorrected + nihtbx_reading_uncorrected  +  lmt_scr_num_correct  +  pea_wiscv_trs
S2 =~ nihtbx_flanker_uncorrected + nihtbx_cardsort_uncorrected  +  nihtbx_pattern_uncorrected
S3 =~ nihtbx_picture_uncorrected  +  pea_ravlt_sd_tc  +  pea_ravlt_ld_tc
G ~~ 0*S1
G ~~ 0*S2
G ~~ 0*S3

S1 ~~ 0*S2
S1 ~~ 0*S3
S2 ~~ 0*S3
nihtbx_reading_uncorrected ~~ 0*nihtbx_reading_uncorrected
'

scores_all = data.frame(subjectkey=data$subjectkey,site_id_l=data$site_id_l)

for (iFold in 1:nFold) {
  print(iFold)
  test_idx = data$site_id_l == u[iFold]
  train_idx = !test_idx
  train_data = data[train_idx,cols]
  test_data = data[test_idx,cols]
  strain = scale(train_data)
  stest = scale(test_data,center=attr(strain,"scaled:center"),scale=attr(strain,"scaled:scale"))
  mdl = cfa(bifactor,data=strain)
  scorestrain = lavPredict(mdl,strain)
  scorestest = lavPredict(mdl,stest)
  scores = data.frame(G=rep(NA,nrow(data)),S1=rep(NA,nrow(data)),S2=rep(NA,nrow(data)),S3=rep(NA,nrow(data)))
  names(scores) = paste0(c("G","S1","S2","S3"),iFold)
  scores[train_idx,] = scorestrain
  scores[test_idx,] = scorestest
  scores_all = cbind(scores_all,scores)
}



#compare to original
sall = scale(data[,cols])
mdl = cfa(bifactor,data=sall,estimator="ML",std.lv=T)

scoresall = lavPredict(mdl,sall)
scores_all$G = scoresall[,1]
scores_all$S1 = scoresall[,2]
scores_all$S2 = scoresall[,3]
scores_all$S3 = scoresall[,3]

scores_all = cbind(scoresall,scores_all)

write.csv(scores_all,"./ABCD_lavaan_gfactor_loso.csv",row.names=F,quote=F,na="NaN")
