function data = load_data(fname)
    

    X = readtable('data.csv',1);
    S = unique(X(:,1));
    
    for s = 1:length(S)
        ix = X(:,1)==S(s);
        for f = 1:length(F)
            data(s).(F{f}) = X(ix,f);
        end    
    
    for s = 1:length(subs)
        ix = D.subject==subs(s);
        data(s).sub = subs(s);
        data(s).intervention = D.agent_interv(ix);
        data(s).trial_start = D.trial_start(ix);
        data(s).choice = D.mine_key.keys(ix);
        data(s).choice_rt = D.mine_key.rt(ix);
        data(s).choice_end = D.choice_end(ix);
        data(s).isi1_start = D.isi1_start(ix);
        data(s).isi1_end = D.isi1_end(ix);
        data(s).feedback_start = D.feedback_start(ix);
        data(s).r = D.feedback(ix);
        data(s).feedback_end = D.feedback_end(ix);
        data(s).isi2_start = D.isi2_start(ix);
        data(s).isi2_end = D.isi2_end(ix);
        data(s).guess_start = D.inter_guess_start(ix);
        data(s).guess_choice = D.intervention_choice_keys.keys(ix);
        data(s).guess_rt = D.intervention_choice_key.rt(ix);
        data(s).guess_end = D.inter_guess_end(ix);
        data(s).iti_start = D.iti_start(ix);
        data(s).iti_end = D.iti_end(ix);
        data(s).iri_start = D.iri_start(ix);
        data(s).iri_end = D.iri_end(ix);
        data(s).cond = D.condition(ix);
        data(s).run = D.run_num(ix);
        data(s).N = length(data(s).c);
        winprob = [D.mine_prob_win_left(ix) D.mine_prob_win_right(ix)];
        for n=1:data(s).N
            [~,k] = max(winprob(n,:));
            if data(s).c(n) == k
                data(s).acc(n) = 1;
            else
                data(s).acc(n) = 0;
            end
        end
        data(s).acc = mean(data(s).acc);
        

    end
    
    ix = [data.acc]>0.6;
    data = data(ix);
