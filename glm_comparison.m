        clear all
        [ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns()
        load_cached_values = false;
        cached_file = fullfile('results', 'glm_comparison.mat');

        EXPT = optCon_expt;
        EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
        goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

        if load_cached_values
            load(cached_file);
        else

            idx = 0;
% 
%             idx = idx + 1;
%             glm(idx).glmodel = 4;
%             glm(idx).name = 'GLM 4';
%             glm(idx).model = 'RPE*psi';
%             glm(idx).pmods = 'RPE*psi';

%{
            idx = idx + 1;
            glm(idx).glmodel = 6;
            glm(idx).name = 'GLM 6';
            glm(idx).model = 'RPE';
            glm(idx).pmods = 'RPE';
%             
            idx = idx + 1;
            glm(idx).glmodel = 5;
            glm(idx).name = 'GLM 5';
            glm(idx).model = 'RPE & psi';
            glm(idx).pmods = 'RPE & psi';
            
            idx = idx + 1;
            glm(idx).glmodel = 11;
            glm(idx).name = 'GLM 11';
            glm(idx).model = 'psi';
            glm(idx).pmods = 'psi';
            %}

            idx = idx + 1;
            glm(idx).glmodel = 55;
            glm(idx).name = 'GLM 55';
            glm(idx).model = 'asym_sticky_rpe';
            glm(idx).pmods = 'RPE';

            idx = idx + 1;
            glm(idx).glmodel = 6;
            glm(idx).name = 'GLM 6';
            glm(idx).model = 'bayesian';
            glm(idx).pmods = 'RPE';

            idx = idx + 1;
            glm(idx).glmodel = 4;
            glm(idx).name = 'GLM 4';
            glm(idx).model = 'bayesian';
            glm(idx).pmods = 'RPE*psi';

            idx = idx + 1;
            glm(idx).glmodel = 5;
            glm(idx).name = 'GLM 5';
            glm(idx).model = 'bayesian';
            glm(idx).pmods = 'RPE,psi';

%             
%             idx = idx + 1;
%             glm(idx).glmodel = 25;
%             glm(idx).name = 'GLM 25';
%             glm(idx).model = 'RPEpsi & RPE & psi';
%             glm(idx).pmods = 'RPEpsi & RPE & psi';



            
%             idx = idx + 1;
%             glm(idx).glmodel = 20;
%             glm(idx).name = 'GLM 20';
%             glm(idx).model = 'RPE_1_minus_psi';
%             glm(idx).pmods = 'RPE_1_minus_psi';

%             idx = idx + 1;
%             glm(idx).glmodel = 21;
%             glm(idx).name = 'GLM 21';
%             glm(idx).model = 'RPE*1-psi & RPE*psi';
%             glm(idx).pmods = 'RPE*1-psi & RPE*psi';
            
            bics = [];
            for i = 1:numel(glm)
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/dorsal_striatum_RL.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/NAC_unnormalized.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/bilateral_caudate_putamen.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/NAC.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/NAC_GLM6_rpe.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/Ca_Pu_GLM4_rpepsi.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt(), glm(i).glmodel, 'masks/Ca_only.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt_nosmooth(), glm(i).glmodel, 'masks/Ca_only.nii', goodSubjs());
                %bic = ccnl_bic(optCon_expt_nosmooth(), glm(i).glmodel, 'masks/NAC.nii', goodSubjs());
                %bic = ccnl_bic(EXPT, glm(i).glmodel, 'masks/sphere_glm44_losses_adversarial_-_wins_adversarial_-_losses_benevolent_+_wins_benevolent_22_30_-18_r=0mm.nii', goodSubjs());
                bic = ccnl_bic(EXPT, glm(i).glmodel, 'masks/Pu.nii', goodSubjs());
                %bic = ccnl_bic(EXPT, glm(i).glmodel, 'masks/NAC.nii', goodSubjs());
                glm(i).bic = bic;
                bics = [bics bic];
            end

            lme = -0.5 * bics;

            [alpha, exp_r, xp, pxp, bor] = bms(lme);

            %save(cached_file);
        end

        % output table
        %
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
