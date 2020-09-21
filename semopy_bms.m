% functional connectivity with SEM

% do BMS with LME generated by semopy_demo.py
% with data from get_betas_for_tetrad.m

disp('feedback_onset');
%load('semopy_S1_feedback_onset_lmes.mat');
load('semopy_Sam_S1_feedback_onset_lmes.mat');

[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor



%{
% fixed effects LR test

uLogL = sum(logliks(:,4)); % bigger (unrestricted) model
rLogL = sum(logliks(:,3)); % smaller (restricted/nested) model
dof = sum(ks(:,4) - ks(:,3)); % degrees of freedom = diff in # params
[h,p,stat,c] = lratiotest(uLogL, rLogL, dof)
disp(p)

% (random effects) LR tests

for s = 1:size(lmes,1)
    uLogL = logliks(s,3); % bigger (unrestricted) model
    rLogL = logliks(s,4); % smaller (restricted/nested) model
    dof = ks(s,3) - ks(s,4); % degrees of freedom = diff in # params
    [h,p,stat,c] = lratiotest(uLogL, rLogL, dof);
    disp(p)
end
%}



disp('trial_onset');
%load('semopy_S1_trial_onset_lmes.mat');
load('semopy_Sam_S1_trial_onset_lmes.mat');
[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor



%{

disp('residuals');
load('semopy_residuals_lmes.mat');
[alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
pxp
bor


% fixed effects LR test

uLogL = sum(logliks(:,4)); % bigger (unrestricted) model
rLogL = sum(logliks(:,3)); % smaller (restricted/nested) model
dof = sum(ks(:,4) - ks(:,3)); % degrees of freedom = diff in # params
[h,p,stat,c] = lratiotest(uLogL, rLogL, dof)
disp(p)
%}
