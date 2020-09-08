% compare functional connectivity
% run after corr_beta_behavior.m

addpath('/users/mtomov13/matlab/');
startup


EXPT = optCon_expt;
%EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
%goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

goodSubjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 25 26 28 29 30 32 33 34]  % correct S3


whiten = true;
filter = true;
%MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
%MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_515voxels_Sphere10.nii', goodSubjs, whiten, filter);
MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, whiten, filter);
%MTG_control_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_contra_ROI_x=-54_y=-26_z=-12_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
MTG_control_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_contra_ROI_x=-54_y=-26_z=-12_515voxels_Sphere10.nii', goodSubjs, whiten, filter);
Pu4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/Put_Sphere4.nii', goodSubjs, whiten, filter);
NAC4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/NAcc_Sphere4.nii', goodSubjs, whiten, filter);

save functional_connectivity.mat


load functional_connectivity.mat

clear conn_mn;
clear conn_mp;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    c = mean(MTG_control_r{j}, 2);
    n = mean(NAC4_r{j}, 2);
    p = mean(Pu4_r{j}, 2);

    conn_mn(j) = corr(m, n);
    conn_mp(j) = corr(m, p);

    conn_cn(j) = corr(c, n);
    conn_cp(j) = corr(c, p);

    conn_pn(j) = corr(p, n);
end

conn_mn = atanh(conn_mn);
conn_mp = atanh(conn_mp);
conn_cn = atanh(conn_cn);
conn_cp = atanh(conn_cp);
conn_pn = atanh(conn_pn);

[h,p,ci,stat] = ttest(conn_mp);
fprintf('MTG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_mn);
fprintf('MTG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_pn);
fprintf('Putamen <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_mp, conn_mn);
fprintf('MTG <-> Putamen vs. MTG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_pn, conn_mp);
fprintf('Putamen <-> VS vs. MTG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_pn, conn_mn);
fprintf('Putamen <-> VS vs. MTG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);


[h,p,ci,stat] = ttest(conn_cp);
fprintf('contralateral MTG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_cn);
fprintf('contralateral MTG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);


[h,p,ci,stat] = ttest(conn_cp, conn_cn);
fprintf('contralateral MTG <-> Putamen vs. contralateral MTG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);



% link functional connectivity to behavior
% run after corr_beta_behavior



[r, p] = corr(conn_mn', acc');
fprintf('MTG (psi ROI) <--> VS (RPE ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);
r1 = r; p1 = p;

[r, p] = corr(conn_mp', acc');
fprintf('MTG (psi ROI) <--> Pu (RPE*psi ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);
r2 = r; p2 = p;

figure;

subplot(1,2,1);
scatter(conn_mn', acc');
xlabel('STS-VS connectivity');
ylabel('accuracy');
text(0.1, 0.8, sprintf('r = %.2f, p = %.2f', r1, p1));
lsline;

subplot(1,2,2);
scatter(conn_mp', acc');
xlabel('STS-Put connectivity');
ylabel('accuracy');
text(0.1, 0.8, sprintf('r = %.2f, p = %.2f', r2, p2));
lsline;

[r, p] = corr(conn_pn', acc');
fprintf('Putamen (RPE*psi ROI) <--> VS (RPE ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);

[r, p] = corr(conn_cn', acc');
fprintf('contralateral MTG (control ROI) <--> NAC (RPE ROI) connectivity does NOT track behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);

[r, p] = corr(conn_cp', acc');
fprintf('contralateral MTG (contral ROI) <--> Pu (RPE*psi ROI) connectivity does NOT track behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);










%{

%MTG_r = ccnl_get_residuals(EXPT, 11, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, false, false);
%MTG_control_r = ccnl_get_residuals(EXPT, 11, '../Momchil/MTG_contra_x=-54_y=-26_z=-12_0voxels_Sphere10.nii', goodSubjs, false, false);
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
conn = atanh(conn);
[r, p] = corr(conn', acc')

clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_control_r{j}, 2);
    n = mean(NAC_r{j}, 2);
    conn(j) = corr(m, n);
end
disp('contralateral MTG (control ROI) <--> NAC (RPE ROI) connectivity does NOT track behavior');
conn = atanh(conn);
[r, p] = corr(conn', acc')

disp('MTG (psi ROI) <--> Ca (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    n = mean(Ca_r{j}, 2);
    conn(j) = corr(m, n);
end
conn = atanh(conn);
[r, p] = corr(conn', acc')

disp('MTG (psi ROI) <--> Pu (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(MTG_r{j}, 2);
    n = mean(Pu_r{j}, 2);
    conn(j) = corr(m, n);
end
conn = atanh(conn);
[r, p] = corr(conn', acc')

%}
