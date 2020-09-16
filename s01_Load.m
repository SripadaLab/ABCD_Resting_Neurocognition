%%
%utility funcions
addpath /home/slab/users/mangstad/repos/Misc_utils/

Exp = '/net/pepper/ABCD/CIFTI/Scripts/rest_neurocognition/';

%%
%Setup files
DataFile = [Exp '/Data/ABCD_rest.csv'];
CorrTemplate = [Exp '/Data/Gordon_Sub_Cere/[Subject].txt'];
NetsFile = [Exp '/Data/gordon_sub_cere_parcels.csv'];

%%
%setup
dat = readtable(DataFile);

N = size(dat,1);
nROI = 418;
P = (nROI*(nROI-1))/2;

netsfile = readtable(NetsFile);
netsfile = netsfile(1:nROI,:);
nets = netsfile.NetworkNumber;

%%
%load connectomes
featuremat = zeros(N,P);

parfor iSubject = 1:N
    Subject = dat.Subject{iSubject};
    fprintf(1,'%d\n',iSubject);
    file = strrep(CorrTemplate,'[Subject]',Subject);
    tmp = load(file);
    tmp = tmp(1:418,1:418);
    tmp = mc_flatten_upper_triangle(tmp);
    featuremat(iSubject,:) = mc_FisherZ(tmp);
end

%%
%setup folds
u = unique(dat.abcd_site_num);
nFold = numel(u);
fold_site_conversion = [[1:nFold]' u];
folds = zeros(size(dat.abcd_site_num));
sitesize = zeros(nFold,1);
for iFold = 1:nFold
    folds(dat.abcd_site_num==u(iFold)) = iFold;
    sitesize(iFold) = sum(folds==iFold);
end

%%
%setup leave one site out variables to predict

losophenoG = [];
losophenoS1 = [];
losophenoS2 = [];
losophenoS3 = [];
for i = 1:nFold
    name = sprintf('G%d',fold_site_conversion(i,2));
    losophenoG = [losophenoG dat.(name)];
    name = sprintf('S1%d',fold_site_conversion(i,2));
    losophenoS1 = [losophenoS1 dat.(name)];
    name = sprintf('S2%d',fold_site_conversion(i,2));
    losophenoS2 = [losophenoS2 dat.(name)];
    name = sprintf('S3%d',fold_site_conversion(i,2));
    losophenoS3 = [losophenoS3 dat.(name)];
end

mainpheno = [dat.G_lavaan dat.S1_lavaan dat.S2_lavaan dat.S3_lavaan];
nihpheno = [dat.nihtbx_cardsort_uncorrected dat.nihtbx_flanker_uncorrected ...
    dat.nihtbx_list_uncorrected dat.nihtbx_pattern_uncorrected dat.nihtbx_picture_uncorrected ...
    dat.nihtbx_picvocab_uncorrected dat.nihtbx_reading_uncorrected dat.pea_ravlt_sd_tc ...
    dat.pea_ravlt_ld_tc dat.pea_wiscv_trs dat.lmt_scr_num_correct];

%%
%get nuisance variables
u = unique(dat.RaceEthnicity);
re = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.RaceEthnicity,u{i});
    re(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        re(idx,i) = NaN;
    end
end
s = sum(re);
i = isnan(s);
re(isnan(re(:,i)),:) = NaN;
re(:,i) = [];
[~,i] = max(nansum(re));
re(:,5) = [];

u = unique(dat.Gender);
gen = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.Gender,u{i});
    gen(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        gen(idx,i) = NaN;
    end
end
s = sum(gen);
i = isnan(s);
gen(isnan(gen(:,i)),:) = NaN;
gen(:,i) = [];
gen(:,2) = [];

u = unique(dat.HighestParentalEducation);
hpe = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HighestParentalEducation,u{i});
    hpe(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        hpe(idx,i) = NaN;
    end
end
s = sum(hpe);
i = isnan(s);
hpe(isnan(hpe(:,i)),:) = NaN;
hpe(:,i) = [];
[~,i] = max(nansum(hpe));
hpe(:,i) = [];

u = unique(dat.HouseholdIncome);
hi = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HouseholdIncome,u{i});
    hi(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        hi(idx,i) = NaN;
    end
end
s = sum(hi);
i = isnan(s);
hi(isnan(hi(:,i)),:) = NaN;
hi(:,i) = [];
[~,i] = max(nansum(hi));
hi(:,i) = [];

u = unique(dat.HouseholdMaritalStatus);
hms = zeros(N,numel(u));
for i = 1:numel(u)
    idx = strcmp(dat.HouseholdMaritalStatus,u{i});
    hms(idx,i) = 1;
    if (strcmp(u{i},'NaN'))
        hms(idx,i) = NaN;
    end
end
s = sum(hms);
i = isnan(s);
hms(isnan(hms(:,i)),:) = NaN;
hms(:,i) = [];
[~,i] = max(nansum(hms));
hms(:,i) = [];

%normal and expanded nuisance
nuisance = [dat.Age dat.Age.^2 gen dat.fd dat.fd.^2 re];
nuisancefull = [dat.Age dat.Age.^2 gen dat.fd dat.fd.^2 re hpe hi hms];


