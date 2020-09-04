% momchil
function latents = get_rational4_latents(subj_original_idx, run)
       load('results_repro32.mat', 'results', 'data', 'bms_results');

       assert(strcmp(func2str(results(1).likfun), 'lik_asym_sticky_rpe'));
       assert(strcmp(func2str(results(2).likfun), 'lik_rational4'));
       assert(length(data) == 32);
       assert(length(results) == 2);
       assert(length(results(1).latents) == 32);
       assert(bms_results.pxp(2) > 0.9);

       s = find([data.sub] == subj_original_idx);
       assert(length(s) == 1);

       if exist('run', 'var')
           w = data(s).run_num == run;
       else
           w = logical(ones(size(data(s).run_num)));
       end

       fn = fieldnames(results(2).latents);
       latents = struct;
       for k=1:numel(fn)
           latents.(fn{k}) = results(2).latents(s).(fn{k})(w,:);
       end
       latents.run_num = data(s).run_num;
end

