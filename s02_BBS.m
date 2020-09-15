%%
%utility funcions
addpath /home/slab/users/mangstad/repos/Misc_utils/

NumComp = 250;
Perms = 10000;

%%
%full sample

results_g_all = mc_bbs(featuremat,mainpheno,nuisance,folds,NumComp,'LOSOPheno',0);
results_g_all_nih = mc_bbs(featuremat,nihpheno,nuisance,folds,NumComp,'LOSOPheno',0,'Scores',results_g_all.Aa);

%%
%extra nuisance
results_g_all_fn = mc_bbs(featuremat,mainpheno,nuisancefull,folds,NumComp,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_nih_fn = mc_bbs(featuremat,nihpheno,nuisancefull,folds,NumComp,'LOSOPheno',0,'Scores',results_g_all.Aa);

%%
%low motion subset
subset1 = dat.fd<0.2;
u = unique(dat.abcd_site_num(subset1));
nFold_sub1 = numel(u);
fold_site_conversion_sub3 = [[1:nFold_sub1]' u];
folds_sub1 = zeros(size(dat.abcd_site_num(subset1)));
sitesize_sub1 = zeros(nFold_sub1,1);
for iFold = 1:nFold_sub1
    folds_sub1(dat.abcd_site_num(subset1)==u(iFold)) = iFold;
    sitesize_sub1(iFold) = sum(folds_sub1==iFold);
end
results_g_lowmotion = mc_bbs(featuremat(subset1,:),mainpheno(subset1,:),nuisance(subset1,:),folds_sub1,NumComp,'LOSOPheno',0);
results_g_lowmotion_nih = mc_bbs(featuremat(subset1,:),nihpheno(subset1,:),nuisance(subset1,:),folds_sub1,NumComp,'LOSOPheno',0,'Scores',results_g_lowmotion.Aa);

%%
%white only
subset2 = strcmp(dat.RaceEthnicity,'White');
u = unique(dat.abcd_site_num(subset2));
nFold_sub2 = numel(u);
fold_site_conversion_sub4 = [[1:nFold_sub2]' u];
folds_sub2 = zeros(size(dat.abcd_site_num(subset2)));
sitesize_sub4 = zeros(nFold_sub2,1);
for iFold = 1:nFold_sub2
    folds_sub2(dat.abcd_site_num(subset2)==u(iFold)) = iFold;
    sitesize_sub4(iFold) = sum(folds_sub2==iFold);
end
results_g_white = mc_bbs(featuremat(subset2,:),mainpheno(subset2,:),nuisance(subset2,:),folds_sub2,NumComp,'LOSOPheno',0);
results_g_white_nih = mc_bbs(featuremat(subset2,:),nihpheno(subset2,:),nuisance(subset2,:),folds_sub2,NumComp,'LOSOPheno',0,'Scores',results_g_white.Aa);

%%
%drop site G
results_g_dsg = mc_bbs(featuremat,losophenoG,nuisance,folds,NumComp,'Scores',results_g_all.Aa,'LOSOPheno',1);
results_g_dss1 = mc_bbs(featuremat,losophenoS1,nuisance,folds,NumComp,'Scores',results_g_all.Aa,'LOSOPheno',1);
results_g_dss2 = mc_bbs(featuremat,losophenoS2,nuisance,folds,NumComp,'Scores',results_g_all.Aa,'LOSOPheno',1);
results_g_dss3 = mc_bbs(featuremat,losophenoS3,nuisance,folds,NumComp,'Scores',results_g_all.Aa,'LOSOPheno',1);


%%
%permutations
results_g_all_perms = mc_bbs_perm(featuremat,mainpheno,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_nih_perms = mc_bbs_perm(featuremat,nihpheno,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);

results_g_all_fn_perms = mc_bbs_perm(featuremat,mainpheno,nuisancefull,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_nih_fn_perms = mc_bbs_perm(featuremat,nihpheno,nuisancefull,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);

results_g_lm_perms = mc_bbs_perm(featuremat(subset1,:),mainpheno(subset1,:),nuisance(subset1,:),folds_sub1,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_lowmotion.Aa);
results_g_lm_nih_perms = mc_bbs_perm(featuremat(subset1,:),nihpheno(subset1,:),nuisance(subset1,:),folds_sub1,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_lowmotion.Aa);

results_g_wo_perms = mc_bbs_perm(featuremat(subset2,:),mainpheno(subset2,:),nuisance(subset2,:),folds_sub2,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_white.Aa);
results_g_wo_nih_perms = mc_bbs_perm(featuremat(subset2,:),nihpheno(subset2,:),nuisance(subset2,:),folds_sub2,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_white.Aa);

results_g_all_perms_dsg = mc_bbs_perm(featuremat,losophenoG,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_perms_dss1 = mc_bbs_perm(featuremat,losophenoS1,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_perms_dss2 = mc_bbs_perm(featuremat,losophenoS2,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);
results_g_all_perms_dss3 = mc_bbs_perm(featuremat,losophenoS3,nuisance,folds,NumComp,Perms,'LOSOPheno',0,'Scores',results_g_all.Aa);

%now calculate proper dropsite performance, since LOSOPheno isn't working
%at the moment for permutations




%permute g pheno and measure dmntpn
inmask = zeros(418,418);
for i = 1:max(nets)
    inmask(nets==i,nets==i) = 1;
end
inmask = inmask==1;

dmntpn = zeros(418,418);
idx = [3,5,8,11,12];
for i = 1:numel(idx)
    for j = 1:numel(idx)
        dmntpn(nets==idx(i),nets==idx(j)) = 1;
        dmntpn(nets==idx(j),nets==idx(i)) = 1;
    end
end
dmntpn = (dmntpn - inmask)>0;
inmask = mc_flatten_upper_triangle(inmask);
dmntpn = mc_flatten_upper_triangle(dmntpn);

mdl_true = fitlm([Aall(good1,1:250) nuisance(good1,:)],dat.G_lavaan(good1));
zcons = zscore(coeffall(:,1:250)*mdl_true.Coefficients.Estimate(2:251));

supra = abs(zcons)>2;
supra = supra';
100*sum(inmask)/numel(inmask)
100*sum(dmntpn)/numel(dmntpn)
100*sum(supra.*inmask)/sum(supra)
100*sum(supra.*dmntpn)/sum(supra)
