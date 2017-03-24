%% PSET 6
% CPNS 34231

%% Problem 1 - Single Neuron Models

% Leaky integrate and fire neuron: tau(dV/dt) = -(V-Vrest)+I(t)

% 1.1) Simulate the neural response for A=3. Create a figure that compares 
% the time course of the input current and the output membrane voltage over
% the first 500 ms. What is the firing rate of the simulated neuron? How many
% spikes does it fire per stimulus cycle?

clear all; close all;

% current parameters
Fc = 40; % Hz
amp = 3; % amplitude

% time parameters
dt = 0.001; % time step (1 ms)
dur = 0.5; % duration (500 ms)

% LIF parameters
Vrest = 0; % resting voltage
Vthresh = 1; % threshold voltage
Vspike = 4; % spiking voltage
Vreset = -1; % reset voltage
tau = 0.015; % membrane time constant (15 ms)

% plot input current and output voltage time courses
[V,I,t,~,firing_rate,spikes_per_cycle] = ...
    LIFNeuron(amp,Fc,dt,dur,Vrest,Vthresh,Vspike,Vreset,tau);

figure;
subplot(2,1,1);
plot(t,V);
xlabel('time (ms)'); ylabel('voltage');
title('Voltage time course');
subplot(2,1,2);
plot(t,I);
xlabel('time (ms)'); ylabel('current');
title('Current time course');

display(firing_rate); display(spikes_per_cycle);

% 1.2) Simulate the neural response for 100 amplitude values between A=0 and A=6.
% Save the spike times that are evoked by the stimuli at each amplitude. 
% Create a spike raster plot of the first 500 ms that illustrates the change 
% in spiking patterns as a function of amplitude. Recreate the Johnson figure
% above, showing how impulses per cycle varies with increasing stimulus amplitude.
% In what way is the model successful in recreating natural behavior? 
% In what ways does it differ? Are any new phenomena predicted by the model? 
% What aspects of the real neuron may not be fully captured by the model?

% generate and sort 100 amplitudes between 0 and 6
low = 0; high = 6;
amps = low + rand(1,100)*(high-low);
amps = sort(amps);

spikes = cell(1,length(amps));
ipc = zeros(1,length(amps)); % impulses per cycle

for i = 1:length(amps)
    [~,~,~,spike_times,~,spikes_per_cycle] = ...
        LIFNeuron(amps(i),Fc,dt,dur,Vrest,Vthresh,Vspike,Vreset,tau);
     spikes{i} = spike_times;
     ipc(i) = spikes_per_cycle;
end

% make a raster plot, of increasing amplitude
figure;
RasterPlot(spikes,'b');
xlabel('time (ms)'); ylabel('trial (increasing amplitude)');
title('Spike raster (increasing amplitude)');

% impulse per cycle versus amplitude plot
figure;
scatter(amps,ipc,'filled');
xlabel('amplitude'); ylabel('impulses per cycle');
title('Impulses per cycle vs amplitude');

% The LIF model is able to recreate the general pattern presented in the
% Johnson figure - flat plateaus where impulses per cycle is constant
% over a range of amplitudes separated by portions where impulses per cycle
% increases approximately linearly with amplitude. Similar to experimental
% results, the plateaus are equally spaced, though they are not separated
% by integer values as in the experimental results. Instead, the plateaus 
% increase by 0.5 in the impulses per cycle, and they do not reach the maximum
% two impulses per cycle that is shown in the actual data. Thus, one new
% phenomenon that may be predicted by this model is that plateaus can occur
% at a constant fractionional number of impulses per cycle. In addition, the Johnson
% figure does not show a plateau at 0 impulses per cycle, whereas the LIF
% model does. Also, the second plateau at 0.5 impulses per cycle spans a
% relatively small range of amplitudes in comparison to the experimental
% results which tend to span a larger range. One might also be able to
% argue that the amplitude ranges between plateaus is smaller than for the
% experimental results (i.e. there is a steeper slope between plateaus).
% The probability of a neuron firing can depend on spiking history (such as
% refractoriness), which is not captured by this LIF model. In addition,
% the model does not incorporate noise which is a property generally found
% in neurons.

%% Problem 2 - Classification

clear all; close all;
load('spikes.mat');

tres = [0.0001:0.0001:0.001 0.002:0.001:0.01 0.02:0.01:0.1]; % temporal resolution (s)
cost = 1./tres; % cost values

% to calculate proportion correct for each fiber and cost value
results = cell(1,length(spikes));
correct = nan(length(spikes), length(cost));

for q = 1:length(cost)
    
    for n = 1:length(spikes)
        % collapse the trials
        trials = {};
        stim = [];
        amp = [];
        t = 1;
        for i = 1:length(spikes{n})
            for j = 1:length(spikes{n}{i})
                trials{t} = spikes{n}{i}{j};
                stim(t) = i;
                amp(t) = j;
                t = t + 1;
             end
        end
        
        % compute a distance matrix for all pairs of trials
        distMat = nan(length(trials));
        for i = 1:length(trials)
            s1 = trials{i};
            for j = 1:length(trials)
                 s2 = trials{j};
                 distMat(i,j) = spkd(s1, s2, cost(q));
            end
        end
        
        % find mean distance to trials in each bandpass
        BP = unique(stim);
        mean_dist = nan(length(trials), length(BP));
        for bp = 1:length(BP)
            idxs = (stim == BP(bp));
            for s = 1:length(trials)
                dist = distMat(s,:);
                dist(s) = nan; % don't consider the spike train itself
                mean_dist(s,bp) = nanmean(dist(idxs));
            end
        end
        
        % assign based on bandpass with minimum average distance to trial
        predict = nan(1, length(trials));
        for s = 1:length(trials)    
            DD = mean_dist(s, :);
            predict(s) = (find(DD == min(DD)));     
        end
        
        % calculate the proportion correct
        result{n}(:,q) = predict;
        correct(n,q) = sum(predict == stim) / length(predict);
        
    end
end

% plot the classification acurracy vs temporal resolution
figure;
h = errorbar(tres*1000, mean(correct,1),std(correct, 0, 1));
set(get(h, 'Parent'), 'XScale', 'log');
xlabel('temporal resolution (ms)'); ylabel('accuracy (%)');
title('Classification performance vs temporal resolution');

% This subset of data illustrates the same general pattern as that shown in
% the 2002 paper.

