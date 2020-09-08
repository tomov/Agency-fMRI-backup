% copy of glm_comparison_voxel
% see https://github.com/tomov/Exploration-Data-Analysis/blob/master/glm_bic_bms_CV.m

clear all;

EXPT = optCon_expt();


%[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns()

EXPT = optCon_expt;
EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

maskfile = '../Momchil/Str_manual.nii';


for s = 1:length(goodSubjs)

    subjs = goodSubjs([1:s-1 s+1:length(goodSubjs)]);

    s
    subjs

    bic6 = ccnl_bic_voxel(EXPT, 6, maskfile, subjs);
    bic4 = ccnl_bic_voxel(EXPT, 4, maskfile, subjs);

    lme6 = -0.5 * bic6;
    lme4 = -0.5 * bic4;

    logBF = lme4 - lme6;

    r = nan(2, size(lme4,2));
    pxp = nan(2, size(lme4,2));
    xp = nan(2, size(lme4,2));
    bor = nan(1, size(lme4,2));
    for i = 1:size(lme4, 2)
        i

        lme = [lme4(:,i) lme6(:,i)];
        [alpha_, exp_r_, xp_, pxp_, bor_, g_] = bms(lme);

        r(:,i) = exp_r_;
        pxp(:,i) = pxp_;
        xp(:,i) = xp_;
        bor(i) = bor_;
    end


    %load('glm_comparison_voxel_CV.mat');


    [mask, Vmask] = ccnl_load_mask(maskfile);

    if s == 1
        VS_sum = zeros(size(mask));
        Put_sum = zeros(size(mask));
    end

    map = zeros(size(mask));

    map(mask) = log(r(1,:) ./ r(2,:)); % log posterior
    %map(mask) = pxp(1,:); % log posterior


    [Pu,V] = ccnl_load_mask('masks/Pu.nii');
    V.fname = '../Momchil/Pu_post.nii'; % !!!!!!!!!!!
    Ca = ccnl_load_mask('masks/Ca.nii');
    Na = ccnl_load_mask('masks/NAC.nii');
    Str = Pu | Ca | Na;

    Str_L = Str;
    Str_R = Str;

    L_xs = 1:round(size(Str,1)/2); % one hemisphere
    R_xs = round(size(Str,1)/2)+1:size(Str,1); % one hemisphere
    Str_R(L_xs,:,:) = 0;
    Str_L(R_xs,:,:) = 0;


    map_L = map;
    map_R = map;

    map_R(~Str_L) = 0;
    map_L(~Str_R) = 0;

    map(~Str) = 0;

    % find min and max in L and R
    [mi,I] = min(map_L, [], 'all', 'linear');
    [x,y,z] = ind2sub(size(map_L), I);
    cor_min_L = [x y z];
    mni_min_L = cor2mni([x y z], V.mat);

    [mi,I] = max(map_L, [], 'all', 'linear');
    [x,y,z] = ind2sub(size(map_L), I);
    cor_max_L = [x y z];
    mni_max_L = cor2mni([x y z], V.mat);

    [mi,I] = min(map_R, [], 'all', 'linear');
    [x,y,z] = ind2sub(size(map_R), I);
    cor_min_R = [x y z];
    mni_min_R = cor2mni([x y z], V.mat);

    [mi,I] = max(map_R, [], 'all', 'linear');
    [x,y,z] = ind2sub(size(map_R), I);
    cor_max_R = [x y z];
    mni_max_R = cor2mni([x y z], V.mat);


    r = 4 / 1.5; % 4 mm, bspmview is +1 compared to ours
    min_L_sphere = create_spherical_mask_helper(mask, cor_min_L(1), cor_min_L(2), cor_min_L(3), r, Vmask);
    max_L_sphere = create_spherical_mask_helper(mask, cor_max_L(1), cor_max_L(2), cor_max_L(3), r, Vmask);
    min_R_sphere = create_spherical_mask_helper(mask, cor_min_R(1), cor_min_R(2), cor_min_R(3), r, Vmask);
    max_R_sphere = create_spherical_mask_helper(mask, cor_max_R(1), cor_max_R(2), cor_max_R(3), r, Vmask);


    V.fname = sprintf('../Momchil/VS_Sphere4_manual_min_CV_%d.nii', goodSubjs(s));
    VS = min_L_sphere | min_R_sphere;
    spm_write_vol(V, VS);

    V.fname = sprintf('../Momchil/Put_Sphere4_manual_max_CV_%d.nii', goodSubjs(s));
    Put = max_L_sphere | max_R_sphere;
    spm_write_vol(V, Put);

    VS_sum = VS_sum + VS;
    Put_sum = Put_sum + Put;
end


V.fname = sprintf('../Momchil/VS_Sphere4_manual_min_CV_sum.nii');
spm_write_vol(V, VS_sum);

V.fname = sprintf('../Momchil/Put_Sphere4_manual_min_CV_sum.nii');
spm_write_vol(V, Put_sum);

