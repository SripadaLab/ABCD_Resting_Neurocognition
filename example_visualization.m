%need to add the Misc_utils repo path for visualization code
addpath /home/slab/users/mangstad/repos/Misc_utils

NetsFile = [Exp '/Data/gordon_sub_cere_parcels.csv'];
nROI = 418;

netsfile = readtable(NetsFile);
netsfile = netsfile(1:nROI,:);
nets = netsfile.NetworkNumber;

cons = load([Exp '/Results/all_15_consensus_vector.txt']);

Tasks = {
    'G'
    'S1'
    'S2'
    'S3'
    'nihtbx_cardsort_uncorrected'
    'nihtbx_flanker_uncorrected'
    'nihtbx_list_uncorrected'
    'nihtbx_pattern_uncorrected'
    'nihtbx_picture_uncorrected'
    'nihtbx_picvocab_uncorrected'
    'nihtbx_reading_uncorrected'
    'pea_ravlt_sd_tc'
    'pea_ravlt_ld_tc'
    'pea_wiscv_trs'
    'lmt_scr_num_correct'
};

%visualize a -1/0/1 z score 2 thresholded consensus
plot_jica_component(cons(:,5)',1,1,2,nets,Tasks{5},[1:16]);

%visualize an -1/0/1 original value 0.01 thresholded component only in a
%subset of networks
plot_jica_component(cons(:,5)',1,0,0.01,nets,Tasks{5},[1 2 3 4 7 10 12]);

%visualize -1/0/1 top 10% of connections 
plot_jica_component(cons(:,5)',1,2,0.9,nets,Tasks{5},[1:16]);

%visualize an unthresholded but weighted component
mc_plot_connectome(cons(:,5)',nets);
