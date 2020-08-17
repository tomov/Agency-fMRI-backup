
%perform a contrast within an ROI mask

EXPT = optCon_expt;
%[allSubjects, subjdirs, goodRuns, goodSubjects, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
EXPT.modeldir = '/ncf/gershman/Lab/Hayley/glmOutput/glms_23_25_accurate';  % so it's the 25 accurate subjects
%[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, 2, 'wins_adversarial - losses_adversarial - wins_benevolent + losses_benevolent', 0.001, '+/-', 0.05, 20, 1);   % get the peak voxel coordinates
%[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, 2, 'wins_adversarial - losses_adversarial', 0.001, '+/-', 0.05, 20, 1);   % get the peak voxel coordinates

acc_subjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 23 26 28 29 30 32 33 34]; 
mask_directory = 'masks/caudate_unnormalized.nii'; %leave this commented if you want to run for real
%mask_directory = cor(1,:); %comment this out if you don't want to do a sanity check

betas_wins_adversarial = ccnl_get_beta(EXPT, 2, 'wins_adversarial', mask_directory, acc_subjs);  % get the betas for each condition
betas_losses_adversarial = ccnl_get_beta(EXPT, 2, 'losses_adversarial', mask_directory, acc_subjs);
betas_wins_benevolent = ccnl_get_beta(EXPT, 2, 'wins_benevolent', mask_directory, acc_subjs);
betas_losses_benevolent = ccnl_get_beta(EXPT, 2, 'losses_benevolent', mask_directory, acc_subjs);

betas_wins_adversarial = mean(betas_wins_adversarial,2);
betas_losses_adversarial = mean(betas_losses_adversarial,2);
betas_wins_benevolent = mean(betas_wins_benevolent,2);
betas_losses_benevolent = mean(betas_losses_benevolent,2);


%Now we have the betas from the peak voxel for each condition.

%First sanity test: wins_adversarial - losses_adversarial shows the same t-stat (12.5328) as the one in the slide / bspmview: 
[h, p, ci, stats] = ttest(betas_wins_adversarial - betas_losses_adversarial)


[h, p, ci, stats] = ttest(betas_wins_benevolent - betas_losses_benevolent)

%However, the 4-way interaction is not significant:
[h, p, ci, stats] = ttest((betas_wins_adversarial - betas_losses_adversarial) - (betas_wins_benevolent - betas_losses_benevolent))







