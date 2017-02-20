%% CPNS 34231
% PSET 3

%% Problem #1: Data Visualization
 
% Plot a raster of the first minute of responses from each neuron in a single 
% plot (each row corresponding to a different neuron).

clear all; close all;
load('retinaData.mat');
spikes = getfield(retinaData, 'spikes');

% isolate the neurons of interest by removing other data
neurons = 1:115;
keep = [1 4 11 15 26 51 80 84 96 105];
remove = find(~ismember(neurons,keep));
spikes(remove) = [];

% plot a raster of the first minute of responses
FirstMinRaster(spikes,'b');

%% Problem #2: Spike-Triggered Average

% A) Use the white noise stimulus and the spikes from neuron #1 to write a 
% function that computes a spatiotemporal receptive field. The stimulus is 
% a 40 by 40 grid in which each checkerboard spot is either on or off. The 
% checkerboard pattern changes at a rate of 30 Hz. Your STRF should range 
% cover lags between 300 ms and -10 ms, with a resolution of 10 ms. Plot 
% your STRF on a 6 by 6 grid of subplots.

close all;
frames = getfield(retinaData, 'stimulusFrames'); % stimulus frames
frame_times = getfield(retinaData, 'stimulusFrameTimes'); % stimulus frames times
lag_times = -10:10:300; % milliseconds
n = length(spikes); % number of neurons

% plot the STRF for different lag times
figure; hold on;
title('Spike-Triggered Average for Neuron #1');
STRFs = zeros(40,40,n,length(lag_times)); % hold the STRFs of all neurons
for i = 1:length(lag_times)
    subplot(6,6,i);
    lag = lag_times(i);
    STRFs(:,:,1,i) = STRF(spikes{1},lag,frames,frame_times);
    title(['STRF: ' num2str(lag) ' ms']);
end

% B) Once you have successfully written a function that can compute the STRF 
% for neuron #1, compute and save the STRF for all 10 neurons.

for i = 2:n
    figure; hold on;
    title(['Spike-Triggered Average for Neuron #' num2str(keep(i))]);
    for j = 1:length(lag_times)
        subplot(6,6,j);
        lag = lag_times(j);
        STRFs(:,:,i,j) = STRF(spikes{i},lag,frames,frame_times);
        title(['STRF: ' num2str(lag) ' ms']);
    end
end

%% Problem #3: Receptive Field Models

% Using lsqcurvefit, fit each of your neuron STRFs using a time-varying 2D 
% Gaussian model. You should fit the Gaussian center, the Gaussian width and 
% 32 different amplitudes (one for each time point).

close all;

% define the time-varying 2D Gaussian model
% fit amplitude, center (x,y), and width (x,y), and offset
Gaussian_model = @(x,xdata) x(1)*exp(-((xdata(:,1) - x(2)).^2/(2*(x(4)^2))...
                            + (xdata(:,2) - x(3)).^2/(2*(x(5)^2)))) + x(6);

amplitudes = zeros(n,length(lag_times)); % hold amplitude time courses for each neuron
lag80_model = zeros(n,6); % hold parameters at lag 80 ms for each neuron
widths = zeros(n*length(lag_times),2); % hold receptive field widths (x and y)
  
[X,Y] = meshgrid(1:40,1:40);
x(:,1) = X(:); % x column
x(:,2) = Y(:); % y column

for i = 1:n
   for j = 1:length(lag_times)
       % reshape the data
       data = STRFs(:,:,i,j);
       data = data(:);
       
       % initial guesses for parameters
       offset = 0.5; % guess the offset
       width = 4.5; % guess the width
       [amp, loc] = max(abs(data - 0.5)); % guess the amplitude by max value
       [R, C] = ind2sub(size(STRFs(:,:,i,j)),loc); % guess center by max location
       x0 = [amp C R width width offset];
       
       % run lsqcurvefit and save certain parameters
       param = lsqcurvefit(Gaussian_model,x0,x,data); % model fitted parameters
       amplitudes(i,j) = param(1); % save all amplitudes
       widths((i-1)*32+j,1) = param(4); widths((i-1)*32+j,2) = param(5); % save all widths
       if lag_times(j) == 80
           lag80_model(i,:) = param; % save all parameters for lag 80 ms
       end
   end
end

% For each neuron, use a series of subplots on a single figure to plot: 
% 1) the empirically determined STRF at time lag 80 ms, 
% 2) the receptive field model at time lag 80 ms, and 
% 3) the amplitude time course from 300 ms to -10 ms. 
% Did the Gaussian successfully fit every neuron?

for i = 1:n
   figure;
   subplot(3,1,1);
   imagesc(STRFs(:,:,i,10)); % empirical STRF at 80 ms lag
   title('Empirical STRF (80 ms lag)');
   subplot(3,1,2); 
   model = Gaussian_model(lag80_model(i,:),x); % receptive field model at 80 ms lag
   model = reshape(model, [40 40]); % reshape into 40x40 matrix
   imagesc(model);
   title('Receptive Field Model STRF (80 ms lag)');
   subplot(3,1,3);
   plot(lag_times,amplitudes(i,:),'b'); % the amplitude time course
   xlabel('lag times (ms)'); ylabel('amplitude');
   title('Amplitude time course: 300 ms to -10 ms');
end

% The Gaussian fit the neurons fairly successfully, well-approximating the
% center, width, and amplitude. However, the features of the average stimulus
% outside of this center are not as well-approximated, highlighted by the
% 6th, 7th, 8th, and 10th neuron whose more complex surrounds are not well
% captured. The 6th neuron in particular is not fitted well by the Gaussian
% model. The model depicts this neuron as a ON neuron but when we look at
% the full STRF of this neuron we see it is actually an OFF neuron. This is
% easier to tell at a 140 ms lag than at the 80 ms lag used. As such, the
% center of the model is also off.

% What is the average receptive field width? What percentage of the neurons 
% in this sample show OFF response properties? What is the average time to 
% peak amplitude for OFF cells? What is the average time to peak amplitude 
% for ON cells?

% compute the average width (x and y)
avg_width = [mean(widths(:,1)) mean(widths(:,2))];
display(avg_width);

% compute the average time to peak for OFF and ON cells
OFF = [1 2 4 5 7 8 9];
ON = [3 10];

average_OFF = 0;
for i = 1:length(OFF)
    % calculate index where OFF neuron's amplitude peaks
    [~,index] = max(amplitudes(OFF(i),:));
    % add the lag time at that index
    average_OFF = average_OFF + lag_times(index);
end
average_OFF = average_OFF/length(OFF); 

average_ON = 0;
for i = 1:length(ON)
    % calculate index where ON neuron's amplitude peaks
    [~,index] = max(amplitudes(ON(i),:));
    % add the lag time at that index
    average_ON = average_ON + lag_times(index);
end
average_ON = average_ON/length(ON);

display(average_OFF); display(average_ON);

% The average receptive field width (x and y), averaged across all neurons and all
% lag times, is shown below. 8 neurons (specifically 1,2,4,5,6,7,8,9), or 
% 80% of neurons in this sample show OFF response properties and 2 show ON 
% response properties. The average time to peak amplitude for OFF and ON cells
% is also shown below.
