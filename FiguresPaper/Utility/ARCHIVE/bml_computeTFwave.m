function out = bml_computeTFwave(cfg)

% % ##################
% 
% % FFT parameters
% nWave = length(wavtime);
% nData = EEG.pnts * EEG.trials; % This line is different from above!!
% nConv = nWave + nData - 1;
% 
% % initialize output time-frequency data
% tf = zeros(length(frex),length(times2save));
% 
% % run convolution
% 
% % now compute the FFT of all trials concatenated
% alldata = reshape( EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:) ,1,[]);
% dataX   = fft( alldata ,nConv );
% 
% % loop over frequencies
% for fi=1:length(frex)
%     
%     % create wavelet and get its FFT
%     % the wavelet doesn't change on each trial...
%     wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
%     waveletX = fft(wavelet,nConv);
%     waveletX = waveletX ./ max(waveletX);
%     
%     % now run convolution in one step
%     as = ifft(waveletX .* dataX);
%     as = as(half_wave+1:end-half_wave);
%     
%     % and reshape back to time X trials
%     as = reshape( as, EEG.pnts, EEG.trials );
%     
%     % compute power and average over trials
%     temppow  = mean( abs(as).^2 ,2);
%     tf(fi,:) = 10*log10( temppow(times2saveidx) / mean(temppow(basetimeidx)) );
% end
% 
% % ##################



% FFT parameters
nWave = length(cfg.wav_time);
half_wave = (length(cfg.wav_time)-1)/2;
%frex = logspace(log10(cfg.min_freq),log10(cfg.max_freq),cfg.num_frex);
frex = cfg.min_freq : cfg.max_freq;
cfg.num_frex = numel(frex);

s = logspace(log10(cfg.range_cycles(1)),log10(cfg.range_cycles(end)),cfg.num_frex) ./ (2*pi*frex);

nChannels = height(cfg.elec_topick);

% initialize output time-frequency data
nTrials = numel(cfg.D.trial);
nPoints = size(cfg.D.trial{1},2);
% 
alldata = cellfun(@(x) x(cfg.elec_topick.Channel,:), cfg.D.trial, 'UniformOutput', false);
alldata_unravel = cell2mat(alldata);

nData = nTrials*nPoints;% EEG.pnts * EEG.trials;

nConv = nWave + nData - 1;

tfdec = zeros(length(frex),nPoints,nTrials,nChannels);


for chani = 1 : height(cfg.elec_topick)
    fprintf("Computed %d/%d Channels \n",chani,height(cfg.elec_topick));
    dataX = fft(alldata_unravel(chani,:), nConv);
    
    % loop over frequencies
    for fi=1:length(frex)    
        fprintf("Computed %d/%d Frequencies \n",fi,length(frex));
        % create wavelet and get its FFT
        % the wavelet doesn't change on each trial...
        wavelet  = exp(2*1i*pi*frex(fi).*cfg.wav_time) .* exp(-cfg.wav_time.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX);
        
        % run convolution
        as = ifft(waveletX .* dataX);
        as = as(half_wave+1:end-half_wave);
        % and reshape back to time X trials
        as = reshape( as, nPoints, nTrials);
        % put power data into big matrix
        tfdec(fi,:,:,chani) = as;
    end
end

out.tfdec = tfdec;
out.elec_topick = cfg.elec_topick;
out.frex = frex;
out.s = s;