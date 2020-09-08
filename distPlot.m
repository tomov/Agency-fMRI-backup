
addpath('/Users/momchil/Dropbox/Research/libs/distributionPlot');

load('semopy_residuals_lmes.mat');

figure;
distributionPlot(bics,'histOpt',2); % histOpt=2 works better for uniform distributions than the default

%{
r = rand(1000,1);
rn = randn(1000,1)*0.38+0.5;
rn2 = [randn(500,1)*0.1+0.27;randn(500,1)*0.1+0.73];
rn2=min(rn2,1);rn2=max(rn2,0);
figure
ah(1)=subplot(2,4,1:2);
boxplot([r,rn,rn2])
ah(2)=subplot(2,4,3:4);
distributionPlot([r,rn,rn2],'histOpt',2); % histOpt=2 works better for uniform distributions than the default
set(ah,'ylim',[-1 2])
%--additional options
data = [randn(100,1);randn(50,1)+4;randn(25,1)+8];
subplot(2,4,5)
distributionPlot(data); % defaults
subplot(2,4,6)
distributionPlot(data,'colormap',copper,'showMM',5,'variableWidth',false) % show density via custom colormap only, show mean/std,
subplot(2,4,7:8)
distributionPlot({data(1:5:end),repmat(data,2,1)},'addSpread',true,'showMM',false,'histOpt',2) %auto-binwidth depends on # of datapoints; for small n, plotting the data is useful

%}
