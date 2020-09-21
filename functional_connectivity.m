% compare functional connectivity
% run after corr_beta_behavior.m

addpath('/users/mtomov13/matlab/');
startup


EXPT = optCon_expt;
%EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
%goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

%goodSubjs = [1 2 4 5 6 8 9 10 11 13 15 16 18 20 21 25 26 28 29 30 32 33 34]  % correct S3
goodSubjs = get_goodSubjs('S1');


whiten = true;
filter = true;
%{
%MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
%MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_515voxels_Sphere10.nii', goodSubjs, whiten, filter);
MTG_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, whiten, filter);
%MTG_control_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_contra_ROI_x=-54_y=-26_z=-12_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
MTG_control_r = ccnl_get_residuals(EXPT, 25, '../Momchil/MTG_contra_ROI_x=-54_y=-26_z=-12_515voxels_Sphere10.nii', goodSubjs, whiten, filter);
Pu4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/Put_Sphere4.nii', goodSubjs, whiten, filter);
NAC4_r = ccnl_get_residuals(EXPT, 25, '../Momchil/NAcc_Sphere4.nii', goodSubjs, whiten, filter);
%}



%{
IFG = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_IFG_ROI_x=42_y=24_z=24_72voxels_Sphere10.nii', goodSubjs, whiten, filter);
Ins = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_AI_ROI_x=36_y=16_z=6_83voxels_Sphere10.nii', goodSubjs, whiten, filter);
VS = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_VS_Sphere4.nii', goodSubjs, whiten, filter);
Put = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_Put_Sphere4.nii', goodSubjs, whiten, filter);
%}

% S1
%{
IFG = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_IFG_ROI_x=42_y=24_z=24_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
Ins = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_AI_ROI_x=36_y=16_z=6_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
VS = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_VS_Sphere4.nii', goodSubjs, whiten, filter);
Put = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_Put_Sphere4.nii', goodSubjs, whiten, filter);
IFG_contra = ccnl_get_residuals(EXPT, 25, '../Momchil/S1_contra_IFG_ROI_x=-42_y=24_z=24_33voxels_Sphere4.nii', goodSubjs, whiten, filter);
 

save S1_functional_connectivity.mat

%}

load S1_functional_connectivity.mat

save final_S1_functional_connectivity.mat

for j = 1:length(goodSubjs)
    ifg = mean(IFG{j}, 2);
    ins = mean(Ins{j}, 2);
    vs = mean(VS{j}, 2);
    put = mean(Put{j}, 2);
    ifgc = mean(IFG_contra{j}, 2);

    conn_ifg_vs(j) = corr(ifg, vs);
    conn_ifg_put(j) = corr(ifg, put);

    conn_ifgc_vs(j) = corr(ifgc, vs);
    conn_ifgc_put(j) = corr(ifgc, put);

    conn_ins_vs(j) = corr(ins, vs);
    conn_ins_put(j) = corr(ins, put);

    conn_put_vs(j) = corr(put, vs);
end


conn_ifg_vs = atanh(conn_ifg_vs);
conn_ifg_put = atanh(conn_ifg_put);
conn_ifgc_vs = atanh(conn_ifg_vs);
conn_ifgc_put = atanh(conn_ifg_put);
conn_ins_vs = atanh(conn_ins_vs);
conn_ins_put = atanh(conn_ins_put);
conn_put_vs = atanh(conn_put_vs);

[h,p,ci,stat] = ttest(conn_ifg_put);
fprintf('IFG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_ifg_vs);
fprintf('IFG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_ifgc_put);
fprintf('contra IFG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_ifgc_vs);
fprintf('contra IFG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_put_vs);
fprintf('Putamen <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_ifg_put, conn_ifg_vs);
fprintf('IFG <-> Putamen vs. IFG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);
p_ifg = p;

[h,p,ci,stat] = ttest(conn_put_vs, conn_ifg_put);
fprintf('Putamen <-> VS vs. IFG <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_put_vs, conn_ifg_vs);
fprintf('Putamen <-> VS vs. IFG <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);


[h,p,ci,stat] = ttest(conn_ins_put);
fprintf('Ins <-> Putamen: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);

[h,p,ci,stat] = ttest(conn_ins_vs);
fprintf('Ins <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);


[h,p,ci,stat] = ttest(conn_ins_put, conn_ins_vs);
fprintf('Ins <-> Putamen vs. Ins <-> VS: t(%d) = %.3f, p = %.10f\n', stat.df, stat.tstat, p);
p_ins = p;




figure;
subplot(1,2,1);
m = [mean(conn_ifg_vs) mean(conn_ifg_put)];
se = [std(conn_ifg_vs) std(conn_ifg_put)] / sqrt(length(conn_ifg_vs));
bar(m);
hold on;
errorbar(m, se, 'Linestyle', 'none', 'color', 'black', 'linewidth', 1);
for i = 1:length(conn_ifg_vs)
    plot([conn_ifg_vs(i) conn_ifg_put(i)], 'color', [0.7 0.7 0.7]);
end
y = max(m + se + 0.06);
line([1 2], [y y], 'color', 'black');
text(1.30, y + 0.02, significance(p_ifg));
xticklabels({'VS', 'Put'});
ylabel('connectivity (atanh(r))');
title('IFG seed');

subplot(1,2,2);
m = [mean(conn_ins_vs) mean(conn_ins_put)];
se = [std(conn_ins_vs) std(conn_ins_put)] / sqrt(length(conn_ins_vs));
bar(m);
hold on;
errorbar(m, se, 'Linestyle', 'none', 'color', 'black', 'linewidth', 1);
for i = 1:length(conn_ins_vs)
    plot([conn_ins_vs(i) conn_ins_put(i)], 'color', [0.7 0.7 0.7]);
end
y = max(m + se + 0.06);
line([1 2], [y y], 'color', 'black');
text(1.30, y + 0.02, significance(p_ins));
xticklabels({'VS', 'Put'});
ylabel('connectivity (atanh(r))');
title('Ins seed');



% link functional connectivity to behavior
% run after corr_beta_behavior



[r, p] = corr(conn_ifg_vs', acc');
fprintf('IFG (psi ROI) <--> VS (RPE ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);
r1 = r; p1 = p;

[r, p] = corr(conn_ifg_put', acc');
fprintf('IFG (psi ROI) <--> Pu (RPE*psi ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);
r2 = r; p2 = p;

figure;

subplot(1,2,1);
scatter(conn_ifg_vs', acc');
xlabel('IFG-VS connectivity');
ylabel('accuracy');
text(0.1, 0.8, sprintf('r = %.2f, p = %.2f', r1, p1));
lsline;

subplot(1,2,2);
scatter(conn_ifg_put', acc');
xlabel('IFG-Put connectivity');
ylabel('accuracy');
text(0.1, 0.8, sprintf('r = %.2f, p = %.2f', r2, p2));
lsline;

[r, p] = corr(conn_put_vs', acc');
fprintf('Putamen (RPE*psi ROI) <--> VS (RPE ROI) connectivity tracks behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);

[r, p] = corr(conn_ins_vs', acc');
fprintf('Ins (psi ROI) <--> NAC (RPE ROI) connectivity does NOT track behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);

[r, p] = corr(conn_ins_put', acc');
fprintf('Ins (psi ROI) <--> Pu (RPE*psi ROI) connectivity does NOT track behavior: r(%d) = %.3f, p = %.10f\n', length(acc), r, p);










%{

%IFG_r = ccnl_get_residuals(EXPT, 11, '../Momchil/IFG_ROI_x=54_y=-26_z=-12_155voxels_Sphere10.nii', goodSubjs, false, false);
%IFG_control_r = ccnl_get_residuals(EXPT, 11, '../Momchil/IFG_contra_x=-54_y=-26_z=-12_0voxels_Sphere10.nii', goodSubjs, false, false);
NAC_r = ccnl_get_residuals(EXPT, 6, 'masks/NAC.nii', goodSubjs, false, false);
Pu_r = ccnl_get_residuals(EXPT, 6, 'masks/Pu.nii', goodSubjs, false, false);
Ca_r = ccnl_get_residuals(EXPT, 6, 'masks/Ca.nii', goodSubjs, false, false);


save scratch.mat

clear conn;
for j = 1:length(goodSubjs)
    m = mean(IFG_r{j}, 2);
    n = mean(NAC_r{j}, 2);
    conn(j) = corr(m, n);
end
disp('IFG (psi ROI) <--> NAC (RPE ROI) connectivity tracks behavior');
conn = atanh(conn);
[r, p] = corr(conn', acc')

clear conn;
for j = 1:length(goodSubjs)
    m = mean(IFG_control_r{j}, 2);
    n = mean(NAC_r{j}, 2);
    conn(j) = corr(m, n);
end
disp('Ins (control ROI) <--> NAC (RPE ROI) connectivity does NOT track behavior');
conn = atanh(conn);
[r, p] = corr(conn', acc')

disp('IFG (psi ROI) <--> Ca (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(IFG_r{j}, 2);
    n = mean(Ca_r{j}, 2);
    conn(j) = corr(m, n);
end
conn = atanh(conn);
[r, p] = corr(conn', acc')

disp('IFG (psi ROI) <--> Pu (RPE ROI) connectivity tracks behavior');
clear conn;
for j = 1:length(goodSubjs)
    m = mean(IFG_r{j}, 2);
    n = mean(Pu_r{j}, 2);
    conn(j) = corr(m, n);
end
conn = atanh(conn);
[r, p] = corr(conn', acc')

%}
