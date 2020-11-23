%%
%Setup output file
OutputFile = [Exp '/Results/ABCD_rest_CIFTI_expressions.csv'];
NumComp = 250;

%%
%run PCA over everyone to get 250 components
[coeffall,~,~,~,expall] = pca(featuremat);

mu = mean(featuremat);
x = bsxfun(@minus,featuremat,mu);
Aall = x*coeffall(:,1:NumComp);

%%
%Save out file with first 250 component expressions for lme model in R
subjectkey = dat.subjectkey;
fid = fopen(OutputFile,'w');

fprintf(fid,'subjectkey,');
for i = 1:NumComp
    temp = sprintf('A%03d',i);
    fprintf(fid,'%s,',temp);
end
fprintf(fid,'\n');
for i = 1:numel(subjectkey)
    fprintf(fid,'%s,',subjectkey{i});
    fprintf(fid,'%15.15f,',Aall(i,:));
    fprintf(fid,'\n');
end
fclose(fid);


%%
%run lme in R
%s04_LME.R

%%
%consensus
betas = load([Exp '/Results/mm_betas.csv']);
cons = coeffall(:,1:NumComp)*betas;

save([Exp '/Results/all_15_consensus_vector.txt'],'cons','-ASCII');

%%
%calculate within network and DMN/TPN connections
zcons = zscore(cons);

close all
plot_jica_component(cons(:,1)',1,1,2,nets,'Mixed Model Consensus',[1:16]);
print([Exp '/Results/consensus.pdf'],'-dpdf','-r600');

tmp = zcons(:,1)';
tmp = mc_unflatten_upper_triangle(tmp);
tmp = tmp + tmp';
save([Exp '/Results/consensus_z_square.txt'],'tmp','-ASCII');

close all
plot_jica_component(cons(:,1)',1,1,2,nets,'Mixed Model Consensus',[1:16]);
print([Exp '/Results/consensus.svg'],'-dsvg','-r600');


dmntpn = zeros(418,418);
idx = [3,5,7,8,11,12];
for i = 1:numel(idx)
    for j = 1:numel(idx)
        dmntpn(nets==idx(i),nets==idx(j)) = 1;
        dmntpn(nets==idx(j),nets==idx(i)) = 1;
    end
end

dmntpn = mc_flatten_upper_triangle(dmntpn);

plot_jica_component(dmntpn,1,0,0,nets,'',[1:16]);

supra = abs(zcons(:,1))>2;
supra = supra';
100*sum(dmntpn)/numel(dmntpn)
100*sum(supra.*dmntpn)/sum(supra)

