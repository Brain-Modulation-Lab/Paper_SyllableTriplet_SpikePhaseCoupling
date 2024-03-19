function [r,p,p_perm,fh] = corr_permute(x,y,type,nperms, fig_flag);

if size(x,1) == 1
    x = x';
end

if size(y,1) == 1
    y = y';
end
nonan = ~isnan(x) & ~isnan(y);

x = x(nonan);
y = y(nonan);

nObs = numel(x);

[r,p] = corr(x,y,'Type',type);
%all_perms = perms(1:nObs);

r_perm = nan(1,nperms);

for perm_i = 1 : nperms
    r_perm(perm_i) = corr(x(randperm(nObs)),y,'Type',type);
end

p_perm = (sum(abs(r_perm - mean(r_perm)) >= abs(r - mean(r_perm)))+1) / (nperms+1);


if fig_flag
    fh = figure('visible','on');
else
    fh = figure('visible','off');
end
histogram(r_perm,floor(sqrt(nperms)));
hold on
xline(r)
title(sprintf('R:%1.2f, p_{perm}:%1.3f, type %s, n_{perm} = %d',r,p_perm,type, nperms))
box off
xlabel('R')
ylabel(' Observations [#]')