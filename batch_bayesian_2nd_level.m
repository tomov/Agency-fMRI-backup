% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'/net/rcss2/srv/export/gershman/share_root/Lab/scripts/matlab/Hayley/batch_bayesian_2nd_level_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
