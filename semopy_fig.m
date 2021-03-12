% see semopy_bms.m

figure('position', [1001 1111 479 228]);


load('fig_5B.mat');
descs
pxp
bor

subplot(3,2,1);
bar(pxp);
xticklabels({'SEM 1A', 'SEM 1B'});
ylim([-0.01 1.01]);
xlim([0.5 2.5]);
ylabel('PXP');



load('fig_5C.mat');
descs
pxp
bor

subplot(3,1,2);
bar(pxp);
xticklabels({'SEM 1A', 'SEM 2A', 'SEM 3A', 'SEM 4A'});
ylim([-0.01 1.01]);
xlim([0.5 4.5]);
ylabel('PXP');



load('fig_5D.mat');
descs
pxp
bor

subplot(3,1,3);
bar(pxp);
xticklabels({'SEM 1B', 'SEM 2B', 'SEM 3B', 'SEM 4B'});
ylim([-0.01 1.01]);
xlim([0.5 4.5]);
ylabel('PXP');

print(gcf, '-dpdf', 'semopy_fig.pdf');

% all


%{

disp('feedback_onset');
load('semopy_S1_feedback_onset_lmes_all.mat');

[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

fb_pxp = pxp;



disp('trial_onset');
load('semopy_S1_trial_onset_lmes.mat');
[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

tr_pxp = pxp;




figure;
subplot(2,1,1);
bar(fb_pxp);
xticklabels({'SEM 8', 'SEM 9', 'SEM 10', 'SEM 11', 'SEM 12', 'SEM 13'});
ylim([-0.01 1.01]);
xlim([0.5 6.5]);
ylabel('PXP');
title('feedback onset');

subplot(2,1,2);
bar(tr_pxp);
xticklabels({'SEM 8', 'SEM 9', 'SEM 10', 'SEM 11', 'SEM 12', 'SEM 13'});
ylim([-0.01 1.01]);
xlim([0.5 6.5]);
ylabel('PXP');
title('trial onset');









% Str

disp('feedback_onset');
load('semopy_S1_feedback_onset_lmes_VS_Put.mat');

[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

fb_pxp = pxp;



disp('trial_onset');
load('semopy_S1_trial_onset_lmes_VS_Put.mat');
[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

tr_pxp = pxp;




figure;
subplot(2,1,1);
bar(fb_pxp);
xticklabels({'SEM 1', 'SEM 2', 'SEM 3'});
ylim([-0.01 1.01]);
xlim([0.5 3.5]);
ylabel('PXP');
title('feedback onset');

subplot(2,1,2);
bar(tr_pxp);
xticklabels({'SEM 1', 'SEM 2', 'SEM 3'});
ylim([-0.01 1.01]);
xlim([0.5 3.5]);
ylabel('PXP');
title('trial onset');





% IFG Ins

disp('feedback_onset');
load('semopy_S1_feedback_onset_lmes_IFG_Ins.mat');

[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

fb_pxp = pxp;



disp('trial_onset');
load('semopy_S1_trial_onset_lmes_IFG_Ins.mat');
[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor

tr_pxp = pxp;




figure;
subplot(2,1,1);
bar(fb_pxp);
xticklabels({'SEM 4', 'SEM 5', 'SEM 6', 'SEM 7'});
ylim([-0.01 1.01]);
xlim([0.5 4.5]);
ylabel('PXP');
title('feedback onset');

subplot(2,1,2);
bar(tr_pxp);
xticklabels({'SEM 4', 'SEM 5', 'SEM 6', 'SEM 7'});
ylim([-0.01 1.01]);
xlim([0.5 4.5]);
ylabel('PXP');
title('trial onset');

%}

