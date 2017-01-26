%% PSET 2
% CPNS 34231

%% Ideal Observer Analysis

clear all; close all;
load('S1_Ideal_Observer_Analysis.mat');

% access spikes and stimulus from data struct
spikes = getfield(Data, 'spikes');
stimuli = getfield(Data, 'stimuli');

% isolate data of 700 amplitude for first neuron
indexes1_700 = find(stimuli{1}(:,2) == 700);
spikes1_700 = spikes{1}(indexes1_700);
stimuli1_700 = stimuli{1}(indexes1_700,:);

% isolate data of 700 amplitude for second neuron
indexes2_700 = find(stimuli{2}(:,2) == 700);
spikes2_700 = spikes{2}(indexes2_700);
stimuli2_700 = stimuli{2}(indexes2_700,:);

% 1a. Plot the rasterplots for these data, sorted by stimulus direction.

% Note: the data is already sorted by stimulus direction
ndir = 16; % the number of directions
figure;
RasterPlotDir(spikes1_700,ndir,stimuli1_700(:,3));
axis([0 0.45 0 80]);
xlabel('Time (s)'); ylabel('Trial (increasing direction)');
title('Neuron 1 (color represents direction)');

figure;
RasterPlotDir(spikes2_700,ndir,stimuli2_700(:,3));
axis([0 0.45 0 80]);
xlabel('Time (s)'); ylabel('Trial (increasing direction)');
title('Neuron 2 (color represents direction)');

% 1b. Plot the direction tuning curve of each neuron (mean firing rate as a 
% function of direction, along with the standard error of the mean).

direction = 0:22.5:337.5;
mean_rate = zeros(1,length(direction));
dur = 0.45; % duration of trials (s)

% calculate mean firing rate for each direction
for i = 1:ndir
    % find number of trials in that direction
    dir_trials = find(stimuli1_700(:,3) == direction(i));
    ntrials = length(dir_trials);
    % find total number of spikes for these trials
    nspikes = 0;
    for j = dir_trials(1):dir_trials(end)
        nspikes = nspikes + length(spikes1_700{j});
    end
    % compute the mean firing rate
    mean_rate(i) = nspikes/(ntrials*dur);
end

% compute the standard error
std_meanrate = std(mean_rate);
std_error = std_meanrate/(ndir^(1/2));
err = repelem(std_error,ndir);

figure;
errorbar(direction,mean_rate,err);
axis([0 350 -5 45]);
xlabel('Direction (degrees)'); ylabel('Mean Firing Rate (spikes/second)');
title('Neuron 1: Mean Firing Rate as a Function of Direction');

% reinstantiate the mean_rate vector
mean_rate = zeros(1,length(direction));

% calculate mean firing rate for each direction
for i = 1:ndir
    % find number of trials in that direction
    dir_trials = find(stimuli2_700(:,3) == direction(i));
    ntrials = length(dir_trials);
    % find total number of spikes for these trials
    nspikes = 0;
    for j = dir_trials(1):dir_trials(end)
        nspikes = nspikes + length(spikes2_700{j});
    end
    % compute the mean firing rate
    mean_rate(i) = nspikes/(ntrials*dur);
end

% compute the standard error
std_meanrate = std(mean_rate);
std_error = std_meanrate/(ndir^(1/2));
err = repelem(std_error,ndir);

figure;
errorbar(direction,mean_rate,err);
axis([0 350 -5 45]);
xlabel('Direction (degrees)'); ylabel('Mean Firing Rate (spikes/second)');
title('Neuron 2: Mean Firing Rate as a Function of Direction');

% 1c. Using ideal observer analysis, compute the neurometric function for 
% each neuron as a function of angular difference in direction. Use the 
% preferred direction as the ?standard? direction. What is the direction 
% discrimination threshold of the two neurons? Is there something abnormal 
% about one of the two neurons according to this analysis?

% FOR NEURON 1

pref1 = 270; % preferred direction for the first neuron
figure;
NeurometricFunc(pref1,stimuli1_700(:,3),spikes1_700);
axis([0 300 0.5 1]);
xlabel('Angular difference (degrees)'); ylabel('Proportion correct');
title('Neuron 1: Neurometric function as function of angular difference in direction');

% FOR NEURON 2

pref2 = 270; % preferred direction for the second neuron
figure;
NeurometricFunc(pref2,stimuli2_700(:,3),spikes2_700);
axis([0 300 0.5 1]);
xlabel('Angular difference (degrees)'); ylabel('Proportion correct');
title('Neuron 2: Neurometric function as function of angular difference in direction');

% If we consider threshold performance as 75% correct, the direction 
% discrimination threshold for the first neuron is around 30 degrees and the
% direction discrimination threshold for the second neuron is around 20 degrees.
% The first neuron expresses a unimodal distribution in the mean firing
% rate as a function of direction. It is curious, however, that there is a
% slight increase in firing rate at 22.5 degrees, far from the preferred
% direction and with succeeding directions showing very little firing. This
% could be due to error or some other feature of the stimulus. It is
% reflected in the neurometric curve by a slight dip in the proportion
% correct at 247.5 degree difference from the preferred direction. The second
% neuron expresses a bimodal distribution in the mean firing rate as a
% function of direction, with peaks separated by around 180 degrees,
% suggesting another feature (such as orientation) may be at work. This too
% is reflected in the neurometric function by a dip in the proportion
% correct at 202.5 degree difference from the preferred direction.

% 1d. Bars are always moving in the direction perpendicular to their orientation.
% Therefore, bars moving at 0 and 180 degrees are both orientated at 90 degrees. 
% Resort the data by orientation. Then, using ideal observer analysis, compute 
% the neurometric function for each neuron as a function of angular difference 
% in orientation. Use the preferred orientation as a the ?standard? orientation. 
% What is the angular (orientation) threshold of the two neurons? 

% Note: the preferred orientations were gleaned from the mean directional 
% firing rates, finding the maximum summed rate between directions
% differing by 180 degrees and then adding 90 degrees

% FOR NEURON 1

pref1 = 180; % the preferred orientation for the first neuron
orientations = mod(stimuli1_700(:,3),180) + 90; % orientations

figure;
NeurometricFunc(pref1,orientations,spikes1_700);
axis([0 90 0.3 1]);
xlabel('Angular difference (degrees)'); ylabel('Proportion correct');
title('Neuron 1: Neurometric function as function of angular difference in orientation');

% FOR NEURON 2

pref2 = 157.5; % the preferred orientation for the second neuron
orientations = mod(stimuli2_700(:,3),180) + 90; % orientations

figure;
NeurometricFunc(pref2,orientations,spikes2_700);
axis([0 90 0.3 1]);
xlabel('Angular difference (degrees)'); ylabel('Proportion correct');
title('Neuron 2: Neurometric function as function of angular difference in orientation');

% If we consider threshold performance as 75% correct, the angular
% orientation threshold of the second neuron is around 27.5 degrees.
% Combining this orientation analysis with the direction analysis of the
% previous section and the bimodal distribution of the second neuron's mean
% firing rates, it is likely that this particular neuron attends to
% orientation. On the other hand, the performance of the first neuron
% levels off at 0.5, or the value of chance. Combining this orientation
% analysis with the direction analysis of the previous section and the
% unimodal distribution of the first neuron's mean firing rate, it is
% likely that this particular neuron attends only to direction.

%% Choice Probability

clear all; close all;
load('choiceData.mat');

spikes = getfield(choiceData, 'spikes');
LIP_spikes = spikes{1};
PFC_spikes = spikes{2};
choices = getfield(choiceData, 'behavioralReport');
LIP_choices = choices{1};
PFC_choices = choices{2};

% 2a. For each neuron, create two raster plots: one for each eventual behavioral
% choice of category. On each raster, draw three lines to separate the pre-trial, 
% sample, delay and test periods.

% sort spike trains by each neuron and each behavioral choice
LIP_choice1 = LIP_spikes(LIP_choices == 1);
LIP_choice2 = LIP_spikes(LIP_choices == 2);
PFC_choice1 = PFC_spikes(PFC_choices == 1);
PFC_choice2 = PFC_spikes(PFC_choices == 2);

% create a raster plot for each neuron and each behavioral choice
% with lines separating pre-trail, sample, delay and test periods
figure;
RasterPlot(LIP_choice1,'b');
xlabel('Time (s)'); ylabel('Trial');
title('LIP Neuron: Choice 1');
figure;
RasterPlot(LIP_choice2,'g');
xlabel('Time (s)'); ylabel('Trial');
title('LIP Neuron: Choice 2');
figure;
RasterPlot(PFC_choice1,'b');
xlabel('Time (s)'); ylabel('Trial');
title('PFC Neuron: Choice 1');
figure;
RasterPlot(PFC_choice2,'g');
xlabel('Time (s)'); ylabel('Trial');
title('PFC Neuron: Choice 2');

% 2b. Calculate and report choice probability values for each of the four 
% periods, for each neuron.

% FOR THE LIP NEURON
disp('LIP neuron');

% determine which choice is preferred
% (i.e. the choice with the higher mean firing rate)
nchoice1 = length(LIP_choice1);
spikes1 = 0;
for i = 1: nchoice1
    spikes1 = spikes1 + length(LIP_choice1{i});
end
rate1 = spikes1/nchoice1;

nchoice2 = length(LIP_choice2);
spikes2 = 0;
for i = 1: nchoice2
    spikes2 = spikes2 + length(LIP_choice2{i});
end
rate2 = spikes2/nchoice2;

display(rate1); display(rate2);

% Comparing rate1 and rate2 shows that 1 is the preferred choice for the LIP neuron
ChoiceProb(LIP_choice1,LIP_choice2);

% FOR THE PFC NEURON
disp('PFC neuron');

% determine which choice is preferred
% (i.e. the choice with the higher mean firing rate)
nchoice1 = length(PFC_choice1);
spikes1 = 0;
for i = 1: nchoice1
    spikes1 = spikes1 + length(PFC_choice1{i});
end
rate1 = spikes1/nchoice1;

nchoice2 = length(PFC_choice2);
spikes2 = 0;
for i = 1: nchoice2
    spikes2 = spikes2 + length(PFC_choice2{i});
end
rate2 = spikes2/nchoice2;

display(rate1); display(rate2);

% Comparing rate1 and rate2 shows that 2 is the preferred choice for the PFC neuron
ChoiceProb(PFC_choice2,PFC_choice1);

% 2c. What is the choice probability in the pre-trial interval? What does this mean?

% The choice probability in the pre-trial interval for the LIP neuron
% (using 1 as the preferred choice) is 0.63258, and the choice probability
% in the pre-trial interval for the PFC neuron (using 2 as the preferred
% choice) is 0.44605. The pre-trial choice probability for the LIP neuron
% is notably higher than chance, meaning that there is a tendency for the
% firing rate of the preferred choice condition to be higher before
% the stimulus is presented. This could mean that the firing rate in
% the LIP neuron before stimulus presentation impacts the eventual choice,
% with a higher pre-trial firing rate tending toward a choice of 1.
