
EXPT = optCon_expt();


%[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns()

EXPT = optCon_expt;
%EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

maskfile = 'masks/NAC.nii';
maskfile = 'masks/Pu.nii';
maskfile = 'masks/mask.nii';

%{
bic55 = ccnl_bic_voxel(EXPT, 55, maskfile, goodSubjs);
bic6 = ccnl_bic_voxel(EXPT, 6, maskfile, goodSubjs);
bic4 = ccnl_bic_voxel(EXPT, 4, maskfile, goodSubjs);
bic5 = ccnl_bic_voxel(EXPT, 5, maskfile, goodSubjs);
%}
bic61 = ccnl_bic_voxel(EXPT, 61, maskfile, goodSubjs);
bic60 = ccnl_bic_voxel(EXPT, 60, maskfile, goodSubjs);

%{
lme55 = -0.5 * bic55;
lme6 = -0.5 * bic6;
lme4 = -0.5 * bic4;
lme5 = -0.5 * bic5;
%}
lme61 = -0.5 * bic61;
lme60 = -0.5 * bic60;

%logBF = lme4 - lme6;
logBF = lme61 - lme60;
logGBF = sum(logBF,1);

: = nan(2, size(lme61,2));
pxp = nan(2, size(lme61,2));
xp = nan(2, size(lme61,2));
bor = nan(1, size(lme61,2));
for i = 1:size(lme61, 2)
    i

    lme = [lme61(:,i) lme60(:,i)];
    [alpha_, exp_r_, xp_, pxp_, bor_, g_] = bms(lme);

    r(:,i) = exp_r_;
    pxp(:,i) = pxp_;
    xp(:,i) = xp_;
    bor(i) = bor_;
end

%save('glm_comparison_voxel.mat', '-v7.3');


%load('glm_comparison_voxel.mat');

%[h,p,ci,stat] = ttest(bic6, bic55);
%stat

mask = ccnl_load_mask(maskfile);

map = zeros(size(mask));
%map(mask) = stat.tstat;
%map(mask) = logGBF;


map(mask) = log(r(1,:) ./ r(2,:)); % log posterior
%map(mask) = pxp(1,:); % log posterior


[Pu,V] = ccnl_load_mask('masks/Pu.nii');
V.fname = '../Momchil/Pu_post.nii'; % !!!!!!!!!!!
Ca = ccnl_load_mask('masks/Ca.nii');
Na = ccnl_load_mask('masks/NAC.nii');
Str = Pu | Ca | Na;

%{
Pu_y = squeeze(sum(sum(Pu,1),3));
miny = min(find(Pu_y));
maxy = max(find(Pu_y));
midy = round((miny * 0.7 + maxy * 0.3) );
Pu(:,midy:end,:) = 0;

spm_write_vol(V, Pu);
%}

%map(~Str) = 0;

bspmview_wrapper(map);
%bspmview_wrapper(Pu);
