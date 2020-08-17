% For each mask (ROI),
%     computes BIC for each model for each subject
%     populates the log model evidences array (rows = subjects, columns =
%     models) which gets passed to BMS
%     which estimates the excedence probabilities (????)
% 
%
[allSubjects, subjdirs, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
 % which ROIs to look at (one at a time) as paths to nifti files
 roi_masks = {fullfile('masks', 'NAC.nii')}; %%%use this one?
%         {fullfile('masks', 'mask.nii'), ...
%          fullfile('masks', 'hippocampus.nii'), ...

%          fullfile('masks', 'vmpfc.nii'), ...
%          fullfile('masks', 'bg.nii'), ...
%          fullfile('masks', 'pallidum.nii'), ...
%          fullfile('masks', 'v1.nii'), ...
%          fullfile('masks', 'visual.nii'), ...
%          fullfile('masks', 'motor.nii'), ...
%          fullfile('masks', 'sensory.nii')};
        
roi_masknames = {};
for i = 1:numel(roi_masks)
    [~, roi_masknames{i}, ~] = fileparts(roi_masks{i});
end
     
% which subjects to analyze (as indices of the subjects array returned by contextGetSubjectsDirsAndRuns)
%subjects = getGoodSubjects();

subjects = goodSubjs();

% which models to consider (as the glmodel value passed to context_create_multi)

%glmodels = [4 55]; % case number for your glms

glmodels = [6 55]; % case number for your glms


alphas = [];
exp_rs = [];
xps = [];
pxps = [];
bors = [];
bics = [];

for roi=roi_masks
    lme = []; % log model evidence
    row_bics = [];
    for glmodel=glmodels
        bic = ccnl_bic(optCon_expt(), glmodel, roi{1}, subjects); %change to ccnl_bic_wholebrain
        lme = [lme, -0.5 * bic];
        row_bics = [row_bics, bic];
    end
    bics = [bics; row_bics];
        
    [alpha,exp_r,xp,pxp,bor,g] = bms(lme);
    alphas = [alphas; alpha];
    exp_rs = [exp_rs; exp_r];
    xps = [xps; xp];
    pxps = [pxps; pxp];
    bors = [bors; bor];
    g
end

roi_masknames{1} = 'NAC';
XP_T = array2table(xps, 'RowNames', roi_masknames, 'VariableNames', {'RPE_PSI_GLM4', 'RPE_TD_GLM55'})

