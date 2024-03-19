function varargout = get_SPCmap(Pairs,SPCvar, T,Tres,cfg,EvtType, PairType, DBType,nPerms)

if isempty(SPCvar)
    SPCvar = {'PPCz'}; % default
end

if ischar(SPCvar)
    SPCvar = cellstr(SPCvar);
end

nSPCvar = numel(SPCvar);
varargout = cell(1, nSPCvar);

nPairs = numel(Pairs.PPCmap);

PairLabel = PairType;

if ~exist('PairType','var') || strcmpi(PairType,'all')
    pair_list = 1 : nPairs;
elseif strcmpi(PairType,'OnlySign')
    pair_list = find([Pairs.Location.nSign]>0)';
elseif strcmpi(PairType,'Top')
    [~, idxUniquePairs] = getUniquePairs(Pairs.Location);
    pair_list = getTop(Pairs,idxUniquePairs);
else % unit type  ['D','N','I','M'];
    
    switch PairType(1)
        % Firing rate mod
        case 'D'
            Unit_Info = get_FRTypes(Pairs);
            unit_type = -1;
            pair_list_tmp = (Unit_Info == unit_type);
        case 'N'
            Unit_Info = get_FRTypes(Pairs);
            unit_type = 0;
            pair_list_tmp = (Unit_Info == unit_type);
        case 'I'
            Unit_Info = get_FRTypes(Pairs);
            unit_type = 1;
            pair_list_tmp = (Unit_Info == unit_type);
            
        case 'M'
            Unit_Info = get_FRTypes(Pairs);
            unit_type = 2;
            pair_list_tmp = (Unit_Info == unit_type);

        % encoding 
        case 'C'
            Unit_Info = get_EncodingTypes(Pairs.Location,'consonant');
            pair_list_tmp = (Unit_Info == 1);
        case 'V'
            Unit_Info = get_EncodingTypes(Pairs.Location,'vowel');
            pair_list_tmp = (Unit_Info == 1);            
        case 'P'
            Unit_Info = get_EncodingTypes(Pairs.Location,'position');
            pair_list_tmp = (Unit_Info == 1);

            
    end


    if contains(PairType,'all')
        pair_list = find(pair_list_tmp);
    elseif contains(PairType,'OnlySign')
        pair_list = find(pair_list_tmp & [Pairs.Location.nSign]>0);
    else
        pair_list = find(pair_list_tmp);
    end
    if size(pair_list,2) == 1
        pair_list = pair_list';

    end
end



    
    %if isnumeric(PairType) % get particular type
%     pair_list = PairType;

%     Pair_Type



if ~exist('DBType','var')
    DBType = 'default';
end

if ~exist('nPerms','var') || isempty(nPerms) || nPerms == 0
    flag_stat = false;
else
    flag_stat = true;
end


switch DBType
    case 'default'
        % create a struct

        Pairs_ = Pairs.TimeEvts;
        [Pairs_.PPC] = deal(Pairs.PPCmap{:});
        [Pairs_.PPCz] = deal(Pairs.PPCzmap{:});
%         [Pairs_.ES] = deal(Pairs.ESmap{:});
        [Pairs_.freq] = deal(mean(cfg.plv.freqsOfInt,2)'); % row vector
%         [Pairs_.Phase] = deal(Pairs.Phasemap{:});

        % [PLV_struct_all.] = deal(ES_mat_all{:});

        % interp
        yQ =  logspace(log10(cfg.plv.freqsOfInt(1,1)), log10(cfg.plv.freqsOfInt(end,1)), 5*size(cfg.plv.freqsOfInt,1));
        xQ = T(1) : Tres : T(2);

        [YQ,XQ] = ndgrid(yQ,xQ);

        % apply interpolation
        for SPCvar_i = 1 : nSPCvar
            SPCmap = [];
            for ii  = pair_list
                if mod(ii,200) == 0
                    fprintf('Getting SPC map (%s) %d - %d \n', SPCvar{SPCvar_i}, ii, nPairs)
                end%[freqGrid, timeGrid] = ndgrid(PLV_struct_all(ii).freq,PLV_struct_all(ii).time);
                switch EvtType
                    case 'Cue'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Cue.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                    case 'Speech'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Speech.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                end
                SPCmap = cat(3,SPCmap,F(YQ,XQ));
            end
            %SPCmap = cellfun(@(x,y) x(:,1:83) ,SPCmap,'UniformOutput',false);

            % apply (very local) gaussian smoothing
            varargout{SPCvar_i} = struct();
            varargout{SPCvar_i}.map = imgaussfilt(squeeze(mean(SPCmap,3,'omitnan')),[.5 .5]);
            varargout{SPCvar_i}.timeQ = xQ;
            varargout{SPCvar_i}.freqQ = yQ;

            if flag_stat
                if isfile(fullfile(cfg.PATHS.saveMatFiles,['PPCzStat_pairs-',(PairLabel), '_evt-', (EvtType),'_nperm-',num2str(nPerms),'.mat']))
                    load(fullfile(cfg.PATHS.saveMatFiles,['PPCzStat_pairs-',(PairLabel), '_evt-', (EvtType),'_nperm-',num2str(nPerms),'.mat']),'Stat')
                else
                    Stat = run_Cluster2DMaps(SPCmap,xQ,yQ, cfg.Baseline, nPerms);
                    save(fullfile(cfg.PATHS.saveMatFiles,['PPCzStat_pairs-',(PairLabel), '_evt-', (EvtType),'_nperm-',num2str(nPerms)]),'Stat')
                end
                varargout{SPCvar_i}.Stat = Stat;

            end

        end

    case 'default-locked'
        Pairs_ = Pairs.TimeEvts;
        [Pairs_.PPC] = deal(Pairs.PPCmap{:});
        [Pairs_.freq] = deal(mean(cfg.plv.freqsOfInt,2)'); % row vector
        [Pairs_.PPCz] = deal(Pairs.PPCzmap{:});



        % [PLV_struct_all.] = deal(ES_mat_all{:});

        % interp
        yQ =  logspace(log10(cfg.plv.freqsOfInt(1,1)), log10(cfg.plv.freqsOfInt(end,1)), 5*size(cfg.plv.freqsOfInt,1));
        xQ = T(1) : Tres : T(2);

        [YQ,XQ] = ndgrid(yQ,xQ);

        % apply interpolation
        for SPCvar_i = 1 : nSPCvar
            SPCmap = [];
            for ii  = pair_list
                if mod(ii,200) == 0
                    fprintf('Getting SPC map (%s) %d - %d \n', SPCvar{SPCvar_i}, ii, nPairs)
                end%[freqGrid, timeGrid] = ndgrid(PLV_struct_all(ii).freq,PLV_struct_all(ii).time);
                switch EvtType
                    case 'Cue'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Cue.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                    case 'Speech'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Speech.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                end
                SPCmap = cat(3,SPCmap,F(YQ,XQ));
            end
            %SPCmap = cellfun(@(x,y) x(:,1:83) ,SPCmap,'UniformOutput',false);

            % apply (very local) gaussian smoothing
            varargout{SPCvar_i} = struct();
            varargout{SPCvar_i}.map = imgaussfilt(squeeze(mean(SPCmap,3,'omitnan')),[.5 .5]);
            varargout{SPCvar_i}.timeQ = xQ;
            varargout{SPCvar_i}.freqQ = yQ;
        end


    case 'gamma_env'
        Pairs_ = Pairs.TimeEvts;
        [Pairs_.PPC] = deal(Pairs.PPCmap{:});
        [Pairs_.freq] = deal(mean(cfg.plv.freqsOfInt,2)'); % row vector
       % [Pairs_.PPCz] = deal(Pairs.PPCzmap{:});



        % [PLV_struct_all.] = deal(ES_mat_all{:});

        % interp
        yQ =  logspace(log10(cfg.plv.freqsOfInt(1,1)), log10(cfg.plv.freqsOfInt(end,1)), 5*size(cfg.plv.freqsOfInt,1));
        xQ = T(1) : Tres : T(2);

        [YQ,XQ] = ndgrid(yQ,xQ);

        % apply interpolation
        for SPCvar_i = 1 : nSPCvar
            SPCmap = [];
            for ii  = pair_list
                if mod(ii,200) == 0
                    fprintf('Getting SPC map (%s) %d - %d \n', SPCvar{SPCvar_i}, ii, nPairs)
                end%[freqGrid, timeGrid] = ndgrid(PLV_struct_all(ii).freq,PLV_struct_all(ii).time);
                switch EvtType
                    case 'Cue'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Cue.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                    case 'Speech'
                        F =griddedInterpolant({Pairs_(ii).freq',Pairs.TimeEvts(ii).Speech.time'},Pairs_(ii).(SPCvar{SPCvar_i}));
                end
                SPCmap = cat(3,SPCmap,F(YQ,XQ));
            end
            %SPCmap = cellfun(@(x,y) x(:,1:83) ,SPCmap,'UniformOutput',false);

            % apply (very local) gaussian smoothing
            varargout{SPCvar_i} = struct();
            varargout{SPCvar_i}.map = imgaussfilt(squeeze(mean(SPCmap,3,'omitnan')),[.5 .5]);
            varargout{SPCvar_i}.timeQ = xQ;
            varargout{SPCvar_i}.freqQ = yQ;
        end
end
end


function pair_list = getTop(Pairs,idxUniquePairs)
nUniquePairs = max(idxUniquePairs);
pair_list = nan(1,nUniquePairs);
MaxSPCUnits = cellfun(@(x) mean(x(:),'omitnan'), Pairs.PPCzmap); %  there was max before
for unit_i = 1 : nUniquePairs
    tmp = find(idxUniquePairs == unit_i);
    [~,idxMaxSPCUnits] = max(MaxSPCUnits(tmp),[],'omitnan');
    pair_list(unit_i) = tmp(idxMaxSPCUnits);
end

end