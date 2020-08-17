EXPT = optCon_expt();


%[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns()

EXPT = optCon_expt;
EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!


bic55 = ccnl_bic_voxel(EXPT, 55, 'masks/mask.nii', goodSubjs);
bic6 = ccnl_bic_voxel(EXPT, 6, 'masks/mask.nii', goodSubjs);
bic4 = ccnl_bic_voxel(EXPT, 4, 'masks/mask.nii', goodSubjs);
bic5 = ccnl_bic_voxel(EXPT, 5, 'masks/mask.nii', goodSubjs);

lme55 = -0.5 * bic55;
lme6 = -0.5 * bic6;
lme4 = -0.5 * bic4;
lme5 = -0.5 * bic5;

logBF = lme5 - lme55;
logGBF = sum(logBF,1);

%[h,p,ci,stat] = ttest(bic6, bic55);
%stat

mask = ccnl_load_mask('masks/mask.nii');

map = zeros(size(mask));
%map(mask) = stat.tstat;
map(mask) = logGBF;

bspmview_wrapper(map);
