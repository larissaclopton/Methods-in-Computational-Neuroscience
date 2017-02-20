function FirstMinRaster(spikes,color)

% plots a raster of the first minute of responses

% INPUTS:
% spikes - cell array of spike times
% color - the color for the plot
   
    hold on;

    n = length(spikes); % the number of neurons
    for i = 1:n
        trial = spikes{i}(spikes{i} < 60*10000); % isolate the first minute of responses
        len_trial = length(trial);
        for j = 1:len_trial
           spkx = [trial(j),trial(j)];
           spky = [i-1,(i-1)+0.9];
           line(spkx,spky,'color',color,'LineWidth',0.25);
        end
    end
 
    xlabel('Time (s)'); ylabel('Neuron');
    title('First Minute Responses of Retinal Ganglion Cells');
end
