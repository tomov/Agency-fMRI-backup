% get variables for SEM / causal modeling / functional connectivity
%

addpath('/users/mtomov13/matlab/');
startup

[allSubjects, subjdirs, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();


EXPT = optCon_expt;
EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34];


%MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, false, false);
%MTG_control_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_contra_x=-54_y=-26_z=-12_0voxels_Sphere10.nii', goodSubjs, false, false);
%Pu4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/Put_Sphere4.nii', goodSubjs, false, false);
%NAC4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/NAcc_Sphere4.nii', goodSubjs, false, false);

masks = {'../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', '../Momchil/Put_Sphere4.nii', '../Momchil/NAcc_Sphere4.nii'};
names = {'MTG', 'Put', 'VS'};


clear T_fb;
clear T_tr;
clear T_res;
clear T_model;

for s = 1:length(goodSubjs)
    subj = goodSubjs(s);
    subj_original_idx = subj_original_indices(subj);

    % get brain data
    for i = 1:length(masks)
        mask = masks{i};
        T_fb{s}(:,i) = mean(ccnl_get_beta_series(optCon_expt(), 57, subj, 'feedback_onset_run', mask), 2);
        T_tr{s}(:,i) = mean(ccnl_get_beta_series(optCon_expt(), 62, subj, 'trial_onset_run', mask), 2);

        res = ccnl_get_residuals(EXPT, 25, mask, subj, true, true);
        T_res{s}(:,i) = mean(res{1}, 2);
    end

    % convert to tables
    T_fb{s} = array2table(T_fb{s}, 'VariableNames', names);
    T_tr{s} = array2table(T_tr{s}, 'VariableNames', names);
    T_res{s} = array2table(T_res{s}, 'VariableNames', names);

    latents = get_rational4_latents(subj_original_idx);
    vars = [latents.psi latents.rpe latents.psi.*latents.rpe];
    which = ismember(latents.run_num, find(goodRuns{subj}));
    vars = vars(which,:);
    T_model{s} = array2table(vars, 'VariableNames', {'psi', 'RPE', 'RPEpsi'});

    T_fb{s} = [T_fb{s} T_model{s}];
    T_tr{s} = [T_tr{s} T_model{s}];
end

save('get_betas_for_tetrad.mat', '-v7.3');
    

restoredefaultpath; % something in the path (likely SPM) fucks with strfind
rmpath(genpath('/n/sw/helmod/apps/centos7/Core/spm/'));

for s = 1:length(goodSubjs)
    subj = goodSubjs(s);

    writetable(T_fb{s}, sprintf('feedback_onset/feedback_onset_SPMsubj%d.txt', subj), 'Delimiter', '\t');
    writetable(T_tr{s}, sprintf('trial_onset/trial_onset_SPMsubj%d.txt', subj), 'Delimiter', '\t');
    writetable(T_res{s}, sprintf('residuals/ROI_residuals_SPMsubj%d.txt', subj), 'Delimiter', '\t');

end
