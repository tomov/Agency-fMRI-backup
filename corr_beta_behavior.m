
%correlate betas to behavior

%load('model_output.mat', 'data')

clear all;

load('model_output_2.mat', 'data', 'results')
[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();


EXPT = optCon_expt;
EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!


%masks = {'masks/NAC_left.nii', 'masks/NAC_right.nii'};
%masks = {'masks/sphere_glm11_psi_54_-26_-12_r=4mm.nii'}; %create contralateral sphere with same coordinates (flipping x to negative; want to see that it's no longer significant); also try bigger and smaller sphere to see if finding is robust
%masks = {'masks/sphere_glm6_RPE_-30_40_-16_r=4mm.nii'}; %left OFC
%masks = {'masks/sphere_glm6_RPE_-6_30_-16_r=4mm.nii'}; %rectus_L
%masks = {'masks/sphere_glm6_RPE_-6_30_-16_r=4mm.nii', 'masks/sphere_glm6_RPE_16_30_-8_r=4mm.nii', 'masks/sphere_glm6_RPE_-28_-6_-24_r=4mm.nii', 'masks/sphere_glm6_RPE_-14_52_8_r=4mm.nii'}; %rectus_L
%masks = {'masks/sphere_glm6_RPE_16_30_-8_r=4mm.nii'}; %mPFC
%masks = {'masks/sphere_glm6_RPE_-28_-6_-24_r=4mm.nii'}; %L_Hipp
%masks = {'masks/sphere_glm6_RPE_-14_52_8_r=4mm.nii'}; %frontal sup medial L
%masks = {'masks/sphere_glm6_RPE_12_50_8_r=4mm.nii'}; %cingulate ant R
%masks = {'masks/sphere_glm6_RPE_8_60_-14_r=4mm.nii'}; %frontal med orb R
%masks = {'masks/sphere_glm6_RPE_-6_30_-16_r=0mm.nii'}; %peak voxel rectus_L

%masks{1} = '../Momchil/insula_ROI_x=42_y=18_z=8_86voxels_Sphere10.nii'; % R AI from GLM 11 psi (26 subj)
%masks{1} = '../Momchil/OFC_ROI_x=24_y=28_z=-18_90voxels_Sphere10.nii'; % R OFC from GLM 11 psi (26 subj)
%masks{1} = '../Momchil/MTG_ROI_x=52_y=-18_z=-10_39voxels_Sphere10.nii'; % R MTL from GLM 11 psi (26 subj)
masks{1} = '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii'; % R MTL from GLM 11 psi (23 subj)



% exclude those not in behavioral dataset
%goodSubjs = goodSubjs(ismember(subj_original_indices(goodSubjs), [data.sub]));



b = [];

for i = 1:length(masks)
    
    
    clear b;
    clear bic;
    clear B;

    behavioral_subj_indices = [];
    for j = 1:length(goodSubjs) 
        fmri_subj_idx = goodSubjs(j);
        %tmp = ccnl_get_beta(optCon_expt, 50, 'temporal', masks{i}, subj); 
        tmp = ccnl_get_beta(EXPT, 11, 'psi', masks{i}, fmri_subj_idx); 
        %tmp = ccnl_get_tmap(optCon_expt, 11, 'psi', masks{i}, subj); 
        b(j) = mean(tmp(:)); %mean(tmp,2) will give you a vector of 4 betas (for each run) - can average across the 2 conds here (look and see the histogram for the positivity bias in psych sci paper)
        bic(j) = ccnl_bic(EXPT, 11, masks{i}, fmri_subj_idx);

        %data.sub = original_subject_idx(subj) %something like this/make sure this is right
        % find subject in behavioral data
       s = find([data.sub] == subj_original_indices(fmri_subj_idx));
       assert(length(s) == 1);
       behavioral_subj_indices(j) = s;

       assert(data(s).sub ~= 4);
    end

    all_b = b;

    %data2 = data([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19, 22:30]); %n=25
    %data2 = data([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19, 22:28]); %n=23
    %data2 = data([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19, 22:28]); %n=23
    data2 = data(behavioral_subj_indices);
    latents = results(2).latents(behavioral_subj_indices);
    
    acc = [];
    accrun = [];
    for j = 1:length(data2)
        acc(j) = mean(data2(j).acc);
    end


    disp('psi betas track behavior');
    [r,p] = corr(all_b', acc') %correlate betas with overall accuracy


    figure; scatter(all_b, acc);
    xlabel('neural beta');
    ylabel('accuracy');
    lsline;


        
    %to do:
    %can also correlate t stats? (ccnl_get_tmap) - accounts for variability 
    %correlate with the subject-specific winning model BIC
    %more masks



    %s is spm subject ID
    %extract trial-by-trial behavioral readout of your choice (accuracy, BIC (exploration paper), pos/neg feedback, condition)
    %correlate this readout with the betas (keep in mind whether you want to do
    %spearman or pearson)
    % dont forget that s needs to match up with betas and behavioral variable


end
