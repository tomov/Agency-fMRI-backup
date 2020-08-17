%-----------------------------------------------------------------------
% Job saved on 18-May-2020 13:10:17 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
EXPT = optCon_expt();

% also see https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=ind1509&L=SPM&D=0&P=437388

matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(EXPT.modeldir, 'model44_bayesian', 'subj1', 'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.LogEv = 'Yes';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.first = 'Yes';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.second = 'No';
%matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.gcon = struct('name', {}, 'convec', {});

load(matlabbatch{1}.spm.stats.fmri_est.spmmat{1});
wa = contains(SPM.xX.name, 'wins_adversarial');
la = contains(SPM.xX.name, 'losses_adversarial');
wb = contains(SPM.xX.name, 'wins_benevolent');
lb = contains(SPM.xX.name, 'losses_benevolent');
wa = wa / sum(wa);
la = la / sum(la);
wb = wb / sum(wb);
lb = lb / sum(lb);
con = wa - la - wb + lb;

matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.gcon(1).name = '4way';
matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.gcon(1).convec = con;

