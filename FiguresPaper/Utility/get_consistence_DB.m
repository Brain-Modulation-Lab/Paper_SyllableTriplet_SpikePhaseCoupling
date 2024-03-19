% list of useful (local) functions
function out = get_consistence_DB(from, to)
fromPairs = [from.Pairs.Location.subj_id, ...
    from.Pairs.Location.S_channel from.Pairs.Location.E_channel];

toPairs = [to.Pairs.Location.subj_id, ...
    to.Pairs.Location.S_channel to.Pairs.Location.E_channel];

n_fromPairs = size(fromPairs,1);

out = to;
idx_from2to = nan(1, n_fromPairs);
for pair_i = 1 : n_fromPairs
    idx_from2to(pair_i) = find(strcmpi(toPairs(:,1), fromPairs(pair_i,1)) & strcmpi(toPairs(:,2), fromPairs(pair_i,2)) & strcmpi(toPairs(:,3), fromPairs(pair_i,3)));  
end

out.Pairs.Location = to.Pairs.Location(idx_from2to,:);
out.Pairs.PPCmap = to.Pairs.PPCmap(idx_from2to);
out.Pairs.TimeEvts = to.Pairs.TimeEvts(idx_from2to);
out.Pairs.PPCzmap = cellfun(@(x,y,z) (x - y)./z, out.Pairs.PPCmap,from.Pairs.meanPPCperm, from.Pairs.stdPPCperm,'UniformOutput',false);

end
