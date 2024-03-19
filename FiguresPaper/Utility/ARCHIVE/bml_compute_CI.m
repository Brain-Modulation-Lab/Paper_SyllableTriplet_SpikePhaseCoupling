function CI = bml_compute_CI(x)
% this fucntion compute the confidence interval
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;
end