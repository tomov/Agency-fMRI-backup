
function EXPT = optCon_expt_nosmooth(local)

    % same as optCon_expt but for nonsmoothed data
    % call with ccnl_fmri_glm_nosmooth

    EXPT = optCon_expt();

    EXPT.modeldir = [EXPT.modeldir, 'glmOutput_nosmooth']; %create these directories on the cluster
    EXPT.rsadir = [EXPT.modeldir, 'rsaOutput_nosmooth'];