%% generate fake data for example
rng(1234);

% # of subjects
n = 1000;

% # of basis vectors
k = 100;

% # of brain ROIs
nr = 100;

% # of brain features calculated from ROIs
p = (nr*(nr-1))/2;

%generate some fake network IDs for plot
nets = repmat([1 1 1 2 3 4 4 5 6 6],1,10);

%first generate some basis vectors and a correlated variable
mu = repmat(0,k+1,1);
r = corr(rand(k,1),rand(k,k));
sigma = eye(k+1);
sigma(1,2:end) = r;
sigma(2:end,1) = r;

basis = mvnrnd(mu,sigma,n);

y = basis(:,1);
basis(:,1) = [];

%generate some random weights to transform basis functions to simulated
%edges
weights = randn(p,k);
X = basis*weights';

%generate a quick random nuisance variable
nuis = randn(n,1);



%% Now move on to the actual analysis example

%s01_Load.m sets paths and loads up the data so that we have a data
%matrix, a phenotype of interest, and nuisance variables. We have that from
%above other than paths and folds
addpath ../Misc_utils

fold10 = randsample(10,n,1);

%number of components to use for BBS prediction
NumComp = 50;

%now run the BBS
%In our analysis this is handled in s02_BBS.m for a number of phenotypes
%and subsets as well as the permutations.

results = mc_bbs(X,y,nuis,fold10,NumComp);

%output structure from mc_bbs has both within-fold and mean across folds
%data for correlation, MSE, NMSE, and R2cv, as well as 95% confidence
%intervals around the mean correlation across folds. It also saves the
%per-fold PCA data. If we want to run prediction on another phenotype, we
%can use the stored PCA data like this
%results2 = mc_bbs(X,y2,nuis,fold10,NumComp,'Scores',results.Aa);

%Run 100 permutations
%This output will include fold and mean permuted performance, and
%permutation p-values for each phenotype passed in.
results_perm = mc_bbs_perm(X,y,nuis,fold10,NumComp,100,'Scores',results.Aa);

%To calculate the consensus connectome, we do PCA on the whole data sample
%and fit a single regression model, then calculate a weighted sum of
%components and betas. In the ABCD analysis, this step actually does PCA in
%matlab (s03_FullPCA.m) then saves the weights out to a file, which is
%loaded in and a linear mixed effects model is fit in R (s04_LME.R) in
%order to properly handle random effects for site and family. Then betas
%from that model are saved and loaded back into matlab to generate the
%consesus. Here we'll just use a linear model in matlab for simplicity.
ca = pca(X);

mu = mean(X);
x = bsxfun(@minus,X,mu);
aa = x*ca(:,1:NumComp);

mdl = fitlm([aa nuis],y);
b = mdl.Coefficients.Estimate(2:(NumComp+1));

consensus = ca(:,1:NumComp)*b;

plot_jica_component(consensus',1,1,2,nets,'sample consensus',[1:10]);
