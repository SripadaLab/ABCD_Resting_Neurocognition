
NetNames = {
    'Somatomotor Hand'
    'Somatomotor Mouth'
    'Cingulo-opercular'
    'Auditory'
    'Default'
    ''
    'Visual'
    'Fronto-parietal'
    'Salience'
    'Subcortical'
    'Ventral Attention'
    'Dorsal Attention'
    'Cerebellum'
    'Uncertain'
    'Cingulo-Parietal'
    'Retrosplenial Temporal'
    };

ZThresh = 3;

allpheno = [mainpheno nihpheno];

nihcons = zscore(cons);

%calculate per cell suprathrehsold pos/neg per consensus
supra_pos = NaN*zeros(max(nets),max(nets),size(nihcons,2));
supra_neg = NaN*zeros(max(nets),max(nets),size(nihcons,2));
zall = nihcons;
% %thresholded
% zall(abs(zall)<ZThresh) = 0;
% zall(zall>0) = 1;
% zall(zall<0) = -1;
indicies = zeros(max(nets),max(nets),2);
for i = 1:max(nets)
    for j = i:max(nets)
        idxi = nets==i;
        idxj = nets==j;
        indicies(i,j,1) = i;
        indicies(j,i,1) = i;
        indicies(i,j,2) = j;
        indicies(j,i,2) = j;
        
        mask = zeros(418,418);
        mask(idxi,idxj) = 1;
        mask(idxj,idxi) = 1;
        mask = mc_flatten_upper_triangle(mask)==1;
%         %thresholded
%         tmpp = sum(zall(mask,:)>0);
%         tmpn = sum(zall(mask,:)<0);
        %weighted
        tmpp = sum(zall(mask,:).*(zall(mask,:)>0));
        tmpn = sum(zall(mask,:).*(zall(mask,:)<0));
        supra_pos(i,j,:) = tmpp;
        supra_neg(i,j,:) = tmpn;
    end
end

tmp1 = reshape(supra_pos,256,15);
tmp2 = reshape(supra_neg,256,15);

ind = reshape(indicies,256,2);

mask1 = std(tmp1,[],2)>0;
mask2 = std(tmp2,[],2)>0;

va = var([tmp1(mask1,:);tmp2(mask2,:)],[],2);
%plot(sort(va))

tmp = cumsum(sort(va,'descend'))/sum(va);

f = find(tmp>0.5);

[~,i] = sort(va,'descend');
u = i(1:f(1));
m1 = [tmp1(mask1,:);tmp2(mask2,:)];
ind = [ind(mask1,:);ind(mask2,:)];
[u ind(u,:)]
posneg = sign((u<=sum(mask1))-0.5);

output = m1(u,:);
g = output(:,1);
s1 = output(:,2);
s2 = output(:,3);
s3 = output(:,4);
cardsort = output(:,5);
flanker = output(:,6);
list = output(:,7);
pattern = output(:,8);
picture = output(:,9);
picvocab = output(:,10);
reading = output(:,11);
ravlt_sd = output(:,12);
ravlt_ld = output(:,13);
wiscv = output(:,14);
lmt = output(:,15);
network1 = NetNames(ind(u,1));
network2 = NetNames(ind(u,2));
output = table(network1,network2,posneg,g,s1,s2,s3,cardsort,flanker,list,pattern,picture,picvocab,reading,ravlt_sd,ravlt_ld,wiscv,lmt);
writetable(output,[Exp sprintf('/Results/network_%d_variance.csv',f(1))]);
