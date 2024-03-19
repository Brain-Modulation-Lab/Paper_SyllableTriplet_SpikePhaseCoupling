session_id = 1;
unit_id = 1;
spkSampRate = fs;
basetimes = coding_events.ITI_starts - PreprocessedSession.starts(session_id) + 1;
trialtimes = event0 - PreprocessedSession.starts(session_id);
IFR_test = IFR.trial{session_id}(unit_id,:);

D = SpikeBin.trial{session_id}(unit_id,:);


