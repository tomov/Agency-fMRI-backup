% glm_comparison.m but with CV ROIs, from glm_comparison_voxel_CV.m

clear all
[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns()
load_cached_values = false;
cached_file = fullfile('results', 'glm_comparison.mat');

EXPT = optCon_expt;
EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!


mask_templates = {...
    '../Momchil/VS_Sphere4_manual_min_CV_%d.nii', ...
    '../Momchil/Put_Sphere4_manual_max_CV_%d.nii', ...
};


for m = 1:length(mask_templates)

    for s = 1:length(goodSubjs)

        mask = sprintf(mask_templates{m}, goodSubjs(s));

        idx = 0;

        idx = idx + 1;
        glm(idx).glmodel = 4;
        glm(idx).name = 'GLM 4';
        glm(idx).model = 'bayesian';
        glm(idx).pmods = 'RPE*psi';

        idx = idx + 1;
        glm(idx).glmodel = 6;
        glm(idx).name = 'GLM 6';
        glm(idx).model = 'bayesian';
        glm(idx).pmods = 'RPE';

        
        bics = [];
        for i = 1:numel(glm)
            bic = ccnl_bic(EXPT, glm(i).glmodel, mask, goodSubjs(s));

            bics(s,i) = bic;
        end
    end


    lmes = -0.5 * bics;

    [alpha, exp_r, xp, pxp, bor, g] = bms(lmes);

    % output table
    %
    disp(mask);
    disp('GLM & model & pmods & PXP\\');
    for i = 1:numel(glm)
        glm(i).pxp = pxp(i);
        fprintf('%s & %s & %s & %.4f \\\\ \n', ...
            glm(i).name, ...
            glm(i).model, ...
            glm(i).pmods, ...
            glm(i).pxp);
    end

    pxp
    bor
end

