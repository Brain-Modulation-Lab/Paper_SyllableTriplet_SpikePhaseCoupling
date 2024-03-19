function Unit_Info = get_FRTypes(Pairs)

PairLoc = [Pairs.Location.S_MNI_X Pairs.Location.S_MNI_Y Pairs.Location.S_MNI_Z];
Unit_id = [Pairs.Location.S_channel Pairs.Location.E_channel Pairs.Location.subj_id];
Unit_FRmod = Pairs.Location.S_typeFRmod;
Unit_id = Unit_id(all(~isnan(PairLoc),2),:);
Unit_FRmod = Unit_FRmod(all(~isnan(PairLoc),2),:);

Unit_Info =  Unit_FRmod;
% PairLoc = PairLoc(all(~isnan(PairLoc),2),:);