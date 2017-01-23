%% Problem Set 1
% CPNS 34231

%% Problem 1: Data visualization

clear all; close all;

% 1.1 Create a raster plot of the data in mtSpikeTimes.mat.

load('mtSpikeTimes.mat');
figure;
subplot(2,2,1);
RasterPlot(mtSpikeTimes,'k');
xlabel('Time (s)');
ylabel('neuron');
title('Raster plot');
subplot(2,2,2);
RasterPlot(mtSpikeTimes,'k');
xlabel('Time (s)');
ylabel('neuron');
title('Raster plot');

% 1.2 Create 2 peri-stimulus time histograms to illustrate the effect of bin size. 

binsize = 0.003;
subplot(2,2,3);
PSTH(mtSpikeTimes,binsize);
axis([0 0.6 0 100]);
xlabel('Time (s)'); ylabel('Firing rate');
title(['PSTH binsize: ' num2str(binsize)]);

binsize = 0.006;
subplot(2,2,4);
PSTH(mtSpikeTimes,binsize);
axis([0 0.6 0 100]);
xlabel('Time (s)'); ylabel('Firing rate');
title(['PSTH binsize: ' num2str(binsize)]);

% A larger binsize gives smaller instantaneous firing rates in the PSTH,
% and a smaller binsize results in more fluctuation of the instanteous
% firing rate.

% Too large of a binsize makes it hard to represent the time-dependent
% spike rate as too large of a time window is considered. However, too
% small of a binsize can cause great fluctuation and make it hard to
% determine the spike rate.

%% Problem 2: Spike Train Variability

clear all; close all;
load('mtSpikeTimes.mat');

% 2.1 Choose a method to discriminate the first spike from background spontaneous
% firing. Find the mean latency of neural response (when first spike occurs) and
% its standard deviation. Plot a histogram of first spike latencies.

% create array of first spike latency for all trials
ntrials = length(mtSpikeTimes);
first_spike_latency = zeros(1,ntrials);

for i = 1:ntrials
    
    % The PSTH's in 2.2 show a large increase in spiking, or peak, around 0.1
    % seconds, so it is likely that spikes before this point are simply
    % background noise. Thus, we consider the first spike latency to be the
    % first spike after 0.09 seconds (just before 0.1).
    first_spike_index = find(mtSpikeTimes{i} > 0.09,1);
    if ~isempty(first_spike_index)
        first_spike_latency(i) = mtSpikeTimes{i}(first_spike_index);
    else
        first_spike_latency(i) = NaN;
    end
 
end

% filter out all NaN's
latencies = first_spike_latency(~isnan(first_spike_latency));

mean_latency = mean(latencies);
display(mean_latency);
std_latencies = std(latencies);
display(std_latencies);
figure;
histogram(latencies);
xlabel('latency (s)'); ylabel('frequency');
title('Histogram of first spike latencies');

% 2.2 Calculate the average firing rate of the neuron within a 30-ms window
% after the first spike. Find the standard deviation and the coefficient of 
% variation. What is the Fano factor of spike?count?in this 30-ms period?

% create array of firing rate for each trial
firing_rates = zeros(1,ntrials);
window = 0.030; % 30ms

% for each trial, find spikes in window after first spike
for i = 1:ntrials
    latency = first_spike_latency(i);
    if ~isnan(latency)
        spikes_in_window = find(mtSpikeTimes{i} > latency & ...
            mtSpikeTimes{i} <= (latency + window));
        firing_rates(i) = length(spikes_in_window)/window;
    else
        firing_rates(i) = NaN;
    end
end

% filter out all NaN's
firing_rates = firing_rates(~isnan(firing_rates));

disp('For firing rates in 30ms window after first spike');
avg_rate = mean(firing_rates);
display(avg_rate);
std_rate = std(firing_rates);
display(std_rate);
CV = std_rate/avg_rate;
display(CV);

% find fano factor (var/mean) of spike count
spike_counts = firing_rates*window;
Fano_factor = var(spike_counts)/mean(spike_counts);
display(Fano_factor);

% reinstantiate firing rates array and increase window size
firing_rates = zeros(1,ntrials);
window = 0.300; % 300ms

% for each trial, find spikes in window after first spike
% also retrieve the spike times and compute the ISIs (2.3)
ISIs = [];
for i = 1:ntrials
    latency = first_spike_latency(i);
    if ~isnan(latency)
        spikes_in_window = find(mtSpikeTimes{i} > latency & ...
            mtSpikeTimes{i} <= (latency + window));
        num_spikes = length(spikes_in_window);
        
        % record the firing rate
        firing_rates(i) = num_spikes/window;
        
        % record the ISIs
        spike_times = mtSpikeTimes{i}(spikes_in_window);
        trial_ISIs = diff(spike_times);
        
        % append trial ISIs to vector holding all ISIs
        ISIs = [ISIs trial_ISIs];    
    else
        firing_rates(i) = NaN;
    end
end

% filter out all NaN's
firing_rates = firing_rates(~isnan(firing_rates));

disp('For firing rates in 300ms window after first spike');
avg_rate = mean(firing_rates);
display(avg_rate);
std_rate = std(firing_rates);
display(std_rate);
CV = std_rate/avg_rate;
display(CV);

% find fano factor (var/mean) of spike count
spike_counts = firing_rates*window;
Fano_factor = var(spike_counts)/mean(spike_counts);
display(Fano_factor);

% 2.3 Plot a histogram of interspike intervals (ISIs) across all trials for 300ms
% window used above. Also plot a curve for the expected ISI count for a
% Poisson neuron firing at the average rate calculated for this period.
% Comment on how well the two curves match.

% NOTE: ISIs were computed above in 2.2 in tandem 
% with firing rates for the 300ms window

% plot ISI histogram
figure;
histogram(ISIs);
xlabel('ISI (s)'); ylabel('Frequency');
title('Histogram of ISIs');

% normalize and overlay poisson distribution
figure; hold on;
h = hist(ISIs);
h = h/sum(h); % normalize
bar(h); 
time = 0:0.01:0.25;
poisson_neuron = exp(-avg_rate*time); % pdf of a poisson distribution
plot(time,poisson_neuron,'r');
title('ISI distribution versus Poisson distribution');
legend('ISI distribution', 'Poisson distribution');

% The first figure plots the histogram of ISIs, which by inspection does
% not look like it will match well to a Poisson distribution, as expected
% by the Fano factor calculated above that does not equal 1. The second
% figure normalizes the distribution of ISIs and overlays a poisson
% distribution based on the average rate of firing from the 300ms window
% calculated above, confirming that the shapes of the curves do not match well.

%% Problem 3: Generalized Linear Models

clear all; close all;
load('hipp_data.mat');

% Plot the animal?s movement trajectory with its position and the spike 
% times for each neuron overlaid. What can you say about the spatial firing 
% properties of these two neurons?

figure; hold on;

% plot the movement trajectory
plot(xN,yN,'k');
axis([-1.25 1.25 -1.25 1.25]);
xlabel('X coordinate'); ylabel('Y coordinate');
title('Movement trajectory');

% overlay spike times for each neuron
spike1_indexes = find(spikes);
spike2_indexes = find(spikes2);
scatter(xN(spike1_indexes),yN(spike1_indexes),12,'r','filled');
scatter(xN(spike2_indexes),yN(spike2_indexes),12,'b','filled');
legend('movement trajectory','spike neuron','spike2 neuron');

% The neurons overlap somewhat, though not completely, in the spatial
% regions in which they fire. In general terms, the spatial regions of
% firing for the first neuron (red) are shifted slightly to the right from the
% spatial regions of firing for the second neuron (blue). Both respond
% frequently to coordinates that are negative in both X and Y.

% Use glmfit to fit an exponential linear intensity model to the spikes
% data. Determine if the maximum likelihood estimates are significant. 
% Plot the maximum likelihood model intensity as a function of (x,y) 
% position. How well does this capture the spatial firing properties of the neuron?

[b1,dev1,stats1] = glmfit([xN yN],spikes,'poisson');
disp('Stats.p for initial model');
display(stats1.p);

% create a grid on which to visualize the model
figure; hold on;
[X,Y]=meshgrid(-1:.1:1);

% compute lambda for each point on the grid using the GLM model
lambda = exp(b1(1)+b1(2)*X+b1(3)*Y);
lambda(find(X.^2 + Y.^2>1))=nan;

% plot lambda as a function of position on the grid
h_mesh = mesh(X,Y,lambda,'AlphaData',0);
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
xlabel('X coordinate'); ylabel('Y coordinate');
title('Linear terms only');

% The p-values in stats1.p confirm that the parameter estimates are
% significant. The plot of the maximum likelihood model intensity as a 
% function of (x,y) position for this exponential linear model does a poor
% job at capturing the spatial firing properties of the spike neuron, though it
% at least depicts firing for coordinates that are negative in X and Y.

% One way to improve on this exponential linear model is to add quadratic 
% functions of xN and yN. Write down and fit a GLM model containing all 
% terms of quadratic and smaller order. Determine if the maximum likelihood 
% estimates for this new model are significant. Plot the maximum likelihood 
% model intensity as a function of (x,y) position. Does this improve the 
% description of the spatial firing properties of the neuron? Compute the 
% AIC (dev+2*(number of parameters)) values for the exponentiated linear 
% and quadratic models. Which model describes the data better?

[b2,dev2,stats2] = glmfit([xN yN xN.^2 yN.^2],spikes,'poisson');
disp('Stats.p for model with quadratic terms');
display(stats2.p);

% % create a grid on which to visualize the model
figure; hold on;
[X,Y]=meshgrid(-1:.1:1);

% compute lambda for each point on the grid using the GLM model
lambda = exp(b2(1)+b2(2)*X+b2(3)*Y+b2(4)*X.^2+b2(5)*Y.^2);
lambda(find(X.^2 + Y.^2>1))=nan;
 
% plot lambda as a function of position on the grid
h_mesh = mesh(X,Y,lambda,'AlphaData',0);
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
xlabel('X coordinate'); ylabel('Y coordinate');
title('Quadratic terms added');

% compute the AIC for exponentiated linear and quadtratic models
linear_AIC = dev1 + 2*length(b1);
quadratic_AIC = dev2 + 2*length(b2);

display(linear_AIC);
display(quadratic_AIC);

% The p-values in stats2.p confirm that the maximum likelihood estimates are
% significant. The plot of this new maximum likelihood model intensity as
% a function of (x,y) more accurately depicts the spatial firing properties
% of the spike neuron, highlighting a region in the 3rd quadrant close to
% the X and Y axis. In addition, this second model has a lower AIC than the 
% strictly linear model, implying less information is lost and that this 
% new model describes the data better than the previous model. 

% For p = 0 to 20 (history going back p ms), fit a GLM model. Plot the AIC
% as a function of p. How far back do you need to go for the model to fit 
% optimally?

% array to hold AIC as function of p from 1 to 20
AICs = zeros(1,20);

% loop over history of the first neuron
for p = 1:20
    % fit the model
    [b,dev,~] = glmfit([[xN yN xN.^2 yN.^2] spikes_hist(:,1:p)],spikes,'poisson');
    % compute the AIC
    AICs(p) = dev + 2*length(b);
end

% plot AIC as a function of p
figure;
plot(1:20,AICs);
xlabel('p ms back in history'); ylabel('AIC');
title('AIC versus p ms back in history');

% AIC is minimized when p = 8, thus one should go back 8ms in the history
% of the neuron to optimally fit the model.

% Augment the GLM model above to include past spiking history of the other
% recorded neuron. Are any of the parameters associated with these spiking 
% interactions significant? Which of your models (of all you have fit) is 
% most parsimonious?

% the value for p is determined from the optimal 
% fit in the previous question
p = 8;
AICs_2 = zeros(1,20);

% loop over history of the second neuron
for pp = 1:20
    [b,dev,stats] = glmfit([[xN yN xN.^2 yN.^2] spikes_hist(:,1:p)...
    spikes2_hist(:,1:pp)],spikes,'poisson');
    AICs_2(pp) = dev + 2*length(b);
end

% plot AIC as a function of pp
figure;
plot(1:20,AICs_2);
xlabel('p ms back in second neuron history'); ylabel('AIC');
title('AIC versus p ms back in second neuron history');

disp('Stats.p for model with spike history of both neurons');
display(stats.p);

% The last 20 values in stats.p show that none of the parameters associated
% with the spiking interactions are significant. The most parsimonious
% model is the simplest model with the lower number of parameters but
% maintaining adequate explanatory power, which in this case is the model
% incorporating linear and quadratic terms as well as up to 8ms back in the
% spiking history of the neuron. This model has the lowest AIC, which
% depends on a lower deviance and smaller number of parameters. We see that
% incorporating the spiking history of the second neuron only increases the
% AIC, suggesting a less optimal model.

% Assume we didn?t know that these neurons were tuned to position. Fit a 
% GLM model with only history and network interaction components. Are
% the parameters related to the network interactions now significant? 
% Explain why the interaction terms might be significant only when we fail 
% to model the spatial component of the firing activity.

[b,dev,stats] = glmfit([spikes_hist spikes2_hist],spikes,'poisson');
disp('Stats.p for model with only history and network components');
display(stats.p);

% The last 20 values of stats.p show that some parameters related to
% network interactions are significant when the spatial components are not
% included in the model. It was shown in the movement trajectory graph
% that the neurons overlapped somewhat in their preferred spatial areas of 
% firing, so it makes sense that these interaction terms are significant.
% However, when spatial components are included in the model, any relevant
% information from interaction with the second neuron becomes redundant. 
% Overall, the spatial component is more informative than the interaction 
% with the second neuron, rendering these terms insignificant when spatial
% components are included in the model.


