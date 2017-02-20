function [avg_stimulus] = STRF(spikes,lag,frames,frame_times)

% computes a STRF based on the spike-triggered average of at a given lag

% INPUTS:
% spikes - cell array of spike times for one neuron
% lag - how far to go back in the stimulus
% frames - stimulus frames changing at a rate of 30Hz
% frame_times - times for each stimulus frame

    n = length(spikes); % the number of spikes
    cum_stimulus = zeros(40,40); % the cumulative stimulus over all spikes

    for i = 1:n
        % go back in time by the lag
        t = double(spikes(i) - lag*10);
        % find the stimulus frame time closest to t
        [~,frame] = min(abs(frame_times - t));
        % retrieve the stimulus at that frame
        stimulus_frame = frames(:,:,frame);
        % add this stimulus frame to the cum_stimulus matrix
        cum_stimulus = cum_stimulus + stimulus_frame;
    end

    % finally, average the values at each checkerboard spot
    avg_stimulus = cum_stimulus/n;
    % plot the avg_stimulus matrix
    imagesc(avg_stimulus)
    
end