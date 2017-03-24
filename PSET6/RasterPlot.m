function RasterPlot(spikes,color)

% plots a raster for a cell array of spike times

% INPUTS:
% spikes - cell array of spike times
% color - plot color

    hold on; 
    
    num_trials = length(spikes);
    for i = 1:num_trials
        trial = spikes{i};
        len_trial = length(trial);
        for j = 1:len_trial
           spkx = [trial(j),trial(j)];
           spky = [(i-1),(i-1)+0.9];
           line(spkx,spky,'color',color,'LineWidth',1);
        end
    end
    
    axis([0 500 0 length(spikes)]);
end

