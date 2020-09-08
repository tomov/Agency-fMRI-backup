
%{
% correlating RPE and RPEpsi

[ subjdirs, nRuns, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

clear rs;
for s = 1:length(goodSubjs)
    subj = goodSubjs(s);

    subj_original_indx = subj_original_indices(subj);

    latents = get_rational4_latents(subj_original_indx);

    RPE = latents.rpe;
    RPEpsi = latents.rpe .* latents.psi;

    [r,p] = corr(RPE, RPEpsi);

    rs(s) = r;
end

zs = atanh(rs);
[h,p,ci,stat] = ttest(zs)
%}


%{
% granger causality

% run after functional_connectivity.m

offset = [-10:1:10];

for i = 1:length(offset)
    o = offset(i);

    conn = nan(length(goodSubjs),1);
    for j = 1:length(goodSubjs)
        n = mean(NAC4_r{j}, 2);
        p = mean(Pu4_r{j}, 2);

        ix = 1:length(n);
        ix = ix + o;
        valid = ix > 0 & ix <= length(n);
        ix(~valid) = 1;
        p = p(ix);
        p = p(valid);
        n = n(valid);
        p = p(ix);
        conn(j) = corr(n, p);
    end

    conn = atanh(conn);

    r(i) = mean(conn);
    se(i) = std(conn) / sqrt(length(conn));
end


figure;
errorbar(offset, r, se);

%}


%{
% run after corr_beta_behavior

MTG_r = ccnl_get_residuals(EXPT, 11, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, false, false);
MTG_control_r = ccnl_get_residuals(EXPT, 11, '../Momchil/MTG_contra_x=-54_y=-26_z=-12_0voxels_Sphere10.nii', goodSubjs, false, false);
NAC_r = ccnl_get_residuals(EXPT, 6, 'masks/NAC.nii', goodSubjs, false, false);
Pu_r = ccnl_get_residuals(EXPT, 6, 'masks/Pu.nii', goodSubjs, false, false);
Ca_r = ccnl_get_residuals(EXPT, 6, 'masks/Ca.nii', goodSubjs, false, false);

save scratch.mat

clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    n = mean(NAC_r{j}, 2);
    conn(j) = corr(m, n);
end
disp('MTG (psi ROI) <--> NAC (RPE ROI) connectivity tracks behavior');
[r, p] = corr(conn', acc')

clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_control_r{j}, 2);
    n = mean(NAC_r{j}, 2);
    conn(j) = corr(m, n);
end
disp('contralateral MTG (control ROI) <--> NAC (RPE ROI) connectivity does NOT track behavior');
[r, p] = corr(conn', acc')

disp('MTG (psi ROI) <--> Ca (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    n = mean(Ca_r{j}, 2);
    conn(j) = corr(m, n);
end
[r, p] = corr(conn', acc')

disp('MTG (psi ROI) <--> Pu (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    n = mean(Pu_r{j}, 2);
    conn(j) = corr(m, n);
end
[r, p] = corr(conn', acc')

%}
