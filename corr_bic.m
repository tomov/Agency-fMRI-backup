function correlate_bic = corr_bic(EXPT,mask,subjects)

%example: corr_bic(optCon_expt, 'VS.nii', goodSubjs)

    neural_bic1 = ccnl_bic(EXPT,4,mask,subjects); %this computes BIC for a single model, but I need to compare 2 - how? do this twice, once for each model

    neural_bic2 = ccnl_bic(EXPT,20,mask,subjects); %this computes BIC for a single model, but I need to compare 2 - how? do this twice, once for each model

    neural_bic = [neural_bic1 neural_bic2];
    
    [alpha,exp_r,xp,pxp,bor,neural_g] = bms(-0.5*neural_bic);
    
    neural_bf = log(neural_g(:,1)./(1 - neural_g(:,1)));


    %[bms_results.alpha, bms_results.exp_r, bms_results.xp, bms_results.pxp, bms_results.bor, bms_results.g] = mfit_bms(results,use_bic)
    
    behav_g = bms_results.g;
    
    assert sum(behav_g) == 1;
    
    behav_bf = log(behav_g(:,1)./(1 - behav_g(:,1)));


    [r,p] = corr(neural_bf, behav_bf);