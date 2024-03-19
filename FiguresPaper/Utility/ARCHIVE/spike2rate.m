function rate=spike2rate(cfg)
% ###################################
% Author:  M. Vissani, 2018
% spike_count,dt,time_window,time_slide
% ###################################
    slide_steps = floor(cfg.time_slide/cfg.dt);
    length_steps = floor(length(cfg.spike_count)/slide_steps);
    windows_steps= floor(floor(cfg.time_window/cfg.dt)/slide_steps);

    spike_hist_tmp = cfg.spike_count(1:(length_steps*slide_steps));
    spike_hist_tmp = reshape(spike_hist_tmp,slide_steps,length_steps);
    spike_hist = sum(spike_hist_tmp);
    kernel = zeros(1,length_steps);
    kernel(1:windows_steps) = 1;
    rate_tmp = conv(spike_hist,kernel);
    rate = rate_tmp(windows_steps:length_steps)/cfg.time_window;
end