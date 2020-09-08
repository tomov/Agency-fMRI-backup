rmpath(genpath('/n/sw/helmod/apps/centos7/Core/spm/'));

load('get_betas_for_tetrad.mat');

for s = 1:length(goodSubjs)
    subj = goodSubjs(s);

    %{
    writetable(T_fb{s}, sprintf('feedback_onset/feedback_onset_SPMsubj%d.txt', subj), 'Delimiter', '\t');
    writetable(T_tr{s}, sprintf('trial_onset/trial_onset_SPMsubj%d.txt', subj), 'Delimiter', '\t');
    %}
    writetable(T_res{s}, sprintf('residuals/ROI_residuals_SPMsubj%d.txt', subj), 'Delimiter', '\t');

end
