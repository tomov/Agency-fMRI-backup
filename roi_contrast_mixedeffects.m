EXPT = optCon_expt;
%[allSubjects, subjdirs, goodRuns, goodSubjects, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
EXPT.modeldir = '/ncf/gershman/Lab/Hayley/glmOutput/glms_23_25_accurate';  % so it's the 25 accurate subjects
%[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, 2, 'wins_adversarial - losses_adversarial - wins_benevolent + losses_benevolent', 0.001, '+/-', 0.05, 20, 1);   % get the peak voxel coordinates
%[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = ccnl_extract_clusters(EXPT, 2, 'wins_adversarial - losses_adversarial', 0.001, '+/-', 0.05, 20, 1);   % get the peak voxel coordinates

acc_subjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 23 26 28 29 30 32 33 34]; 
mask_directory = 'masks/dorsal_striatum_RL.nii'; %leave this commented if you want to run for real
%mask_directory = 'masks/Ca_unnormalized.nii';
%mask_directory = cor(1,:); %comment this out if you don't want to do a sanity check

betas_wins_adversarial = ccnl_get_beta_mixedeffects(EXPT, 2, 'wins_adversarial', mask_directory, acc_subjs);  % get the betas for each condition
betas_losses_adversarial = ccnl_get_beta_mixedeffects(EXPT, 2, 'losses_adversarial', mask_directory, acc_subjs);
betas_wins_benevolent = ccnl_get_beta_mixedeffects(EXPT, 2, 'wins_benevolent', mask_directory, acc_subjs);
betas_losses_benevolent = ccnl_get_beta_mixedeffects(EXPT, 2, 'losses_benevolent', mask_directory, acc_subjs);


b = [];
cond = [];
S = [];

for s = 1:length(acc_subjs)
    subj = acc_subjs(s);
    b = [b;betas_wins_adversarial{s}]; 
    b = [b;betas_losses_adversarial{s}]; 
    b = [b;betas_wins_benevolent{s}]; 
    b = [b;betas_losses_benevolent{s}]; 
    cond = [cond;ones(size(betas_wins_adversarial{s}))];%cond1
    cond = [cond;2*ones(size(betas_losses_adversarial{s}))];%cond2
    cond = [cond;3*ones(size(betas_wins_benevolent{s}))];%cond3
    cond = [cond;4*ones(size(betas_losses_benevolent{s}))];%cond4
    S = [S;s*ones(size(betas_wins_adversarial{s}))];
    S = [S;s*ones(size(betas_losses_adversarial{s}))];
    S = [S;s*ones(size(betas_wins_benevolent{s}))];
    S = [S;s*ones(size(betas_losses_benevolent{s}))];
    
end

cond = categorical(cond);


model_table = table(b,cond,S);


formula = 'b ~ -1 +cond+(-1+cond|S)'

results_lme = fitglme(model_table,formula,'Distribution','Normal','Link','Identity','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal', 'DummyVarCoding','Full')
    


H = [1 -1 -1 1];   % contrast matrix for full conditionxvalence interaction
%H = [1 0 -1 0];   % contrast matrix for full wins_adversarial - wins_benevolent
[p,F,DF1,DF2] = coefTest(results_lme,H)
