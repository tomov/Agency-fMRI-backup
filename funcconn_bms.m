
for dirname = {'S1_feedback_onset', 'S1_trial_onset'}

    lmes = [];

    files = dir(dirname{1});
    for i = 1:length(files)
        filename = fullfile(dirname{1}, files(i).name);

        if endsWith(filename, '.txt')
            filename

            tbl = readtable(filename);

            clear res;

            formula = 'Put ~ RPEpsi';
            res{1} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + Ins';
            res{2} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + IFG';
            res{3} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + Ins + IFG';
            res{4} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + Ins + VS';
            res{5} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + IFG + VS';
            res{6} = fitglme(tbl,formula);

            formula = 'Put ~ RPEpsi + Ins + IFG + VS';
            res{7} = fitglme(tbl,formula);

            bic = cellfun(@(x) x.ModelCriterion.BIC, res);
            lmes = [lmes; - 0.5 * bic];
        end
    end


    disp(dirname{1})
    [alpha,exp_r,xp,pxp,bor,g] = bms(lmes);
    pxp
    bor
  
    %formula = 'C ~ -1 + V + RU + (-1 + V + RU|S)';
    %results_VRU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace',  'EBMethod', 'TrustRegion2D','CovariancePattern','diagonal')

end
