function data_vec = vectorize_groupstats(data, field)

n_fields = numel(field);
n_rows = numel(data);

data_vec = cell(1,n_fields);

for field_i = 1 : n_fields
    if isstruct(data)      
        for row_i = 1 : n_rows
            data_vec{1,field_i} = [data_vec{1,field_i}; data(row_i).(field)];
        end
    elseif istable(data)
        data_vec{1,field_i} = [data.(field)];
    end
end

if n_fields == 1
    data_vec = data_vec{1,1};
end



% STA_pow = [];
% % STA_real = [];
% n_pairs = numel(STA_group);
% for pair_i = 1 : n_pairs
%     STA_pow = [STA_pow; STA_group(pair_i).POW];
%     %     STA_real = [STA_real; STA_group(pair_i).REAL];
%     
% end
% STA_pow = [STA_group.POW];