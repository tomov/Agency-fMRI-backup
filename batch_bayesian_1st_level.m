% List of open inputs
rmpath('/n/sw/helmod/apps/centos7/Core/spm/12.7487-fasrc01/external/fieldtrip/compat/matlablt2016b/');

nrun = 1; % enter the number of runs here
jobfile = {'/ncf/gershman/Lab/scripts/matlab/Hayley/batch_bayesian_1st_level_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
