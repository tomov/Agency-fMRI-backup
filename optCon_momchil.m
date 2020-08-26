function EXPT = optCon_momchil()
    
    % used for momchil's GLMs b/c of file permissions

    EXPT = optCon_expt(false);

    EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
    %goodSubjs = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!

    %EXPT.modeldir = fullfile(EXPT.modeldir, 's3_analyses_aug2020');
    %EXPT.modeldir = '/ncf/gershman/Lab/Momchil-Agency/glmOutput';

