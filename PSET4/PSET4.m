%% PSET4
% CPNS 34231

% The file mtNeuron.mat contains the responses of a single directionally tuned
% MT neuron to random dot stimuli moving coherently in directions varying from
% -90 to 90 degrees relative to the previously measured preferred direction of 
% the neuron. The thirteen stimuli were each presented 184 times. Each stimulus 
% began at time 0 and continued for 256 ms. Recordings continued until 512 ms 
% after the beginning of the stimulus. The array has the dimensions 256?13?184. 
% The first dimension is time in 2 ms bins, the second dimension motion direction,
% the third dimension the repeated presentations.

%% Problem 1

% Create a raster of all the spike trains in one plot, sorted by stimulus 
% direction. Rasters corresponding to different directions should be 
% represented in different colors.

clear all; close all;
load('mtNeuron.mat');

spikes = getfield(mtNeuron,'data');
times = getfield(mtNeuron,'time');
figure;
RasterPlotDir(spikes,times,13,184);
xlabel('Time (ms)'); ylabel('Trial (color represents direction)');
title('Raster Plot (direction -90 to 90 degrees)');

%% Problem 2

% Compute and plot the mutual information between cumulative spike count
% and motion direction as a function of time. 

% parameters
[dur,ndir,ntrials] = size(spikes);
p_dir = 1.0/ndir; % probability of a given direction
max_spikes = dur; % the maximum possible number of spikes
cum_spikes = cumsum(spikes,1); % the cumulative spikes

% variables
mutual_information = zeros(1,dur);
pt_n_dir = zeros(dur,ndir,max_spikes);
pt_n = zeros(dur,max_spikes);

for t = 1:dur
    
    % first, compute pt_n_dir for all n and all dir
    for n = 1:max_spikes
        for dir = 1:ndir
            pt_n_dir(t,dir,n) = length(find(cum_spikes(t,dir,:) == (n-1)))/ntrials;
        end
    end
    
    % next, compute pt_n for all n
    for n = 1:max_spikes
        pt_n(t,n) =  sum(p_dir*pt_n_dir(t,:,n));
    end
    
    % use these values to compute the mutual information
    total_result = 0;
    for dir = 1:ndir
        direction_result = 0;
        for n = 1:max_spikes
            if(pt_n(t,n)*pt_n_dir(t,dir,n) ~= 0)
                direction_result = direction_result + pt_n_dir(t,dir,n)*log2(pt_n_dir(t,dir,n)/pt_n(t,n));
            end
        end
        total_result = total_result + p_dir*direction_result;
    end
    mutual_information(t) = total_result;
    
end

figure;
plot(times*1000,mutual_information);
xlabel('Time (ms)'); ylabel('Mutual Information (bits)');
title('Mutual Information between Cumulative Spike Count and Motion Direction');

% Note: This method for computing the mutual information is outlined in
% Osborne et al, 2004.

%% Problem 3

% Determine the latency of this neuron in a principled way. What proportion
% of the mutual information is available within the first 50 ms of the neural
% response? Within the first 100 ms?

% determine the latency of the neuron (in ms)
latency = 0;
counts = zeros(1,dur); % total spike count per time point
for t = 1:dur
    for dir = 1:ndir
        counts(t) = counts(t) + sum(spikes(t,dir,:));
    end
    if t ~= 1
        previous_mean = mean(counts(1:t-1));
        if counts(t) - previous_mean > 6*std(counts(1:t-1))
            latency = 1000*times(t);
            break
        end
    end
end
display(latency);

% compute the proportion of mutual information available after 50 ms and 100 ms
binsize = 0.002; % bin size of the recordings
max_info = max(mutual_information); % the maximum mutual information
proportion50 = mutual_information(((latency+50)/1000)/binsize)/max_info;
proportion100 = mutual_information(((latency+100)/1000)/binsize)/max_info;
display(proportion50); display(proportion100);

% Based on the raster plot, the latency of the neuron appears to be around
% 90 milliseconds as this tends to be where firing picks up significantly for
% directions focused around the preferred direction. A principled way to
% determine this latency is to find the time point where the difference
% between the current spike count (summed across all directions and all
% trials) and the mean spike count of all previous time points has a value
% greater than 6 times the standard deviation of the spike counts from the
% previous time points. This resulted in a latency of 94 milliseconds. 

% Based on this latency, a greater proportion of the mutual information is 
% available within 100 ms of the response than within 50 ms of the response.
% The exact proportions are shown below. This means more "information" can
% be obtained about the neuron from the spike counts within 100 ms of the
% response.

