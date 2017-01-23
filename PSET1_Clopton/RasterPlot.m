function RasterPlot(spikes,color)

% plots a raster for a cell array of spike times

% INPUTS:
% spikes - cell array of spike times
% color - plot color
% fs - fontsize

    hold on; 
    
    num_trials = length(spikes);
    for i = 1:num_trials
        trial = spikes{i};
        len_trial = length(trial);
        for j = 1:len_trial
           spkx = [trial(j),trial(j)];
           spky = [i,i+0.9];
           line(spkx,spky,'color',color,'LineWidth',1);
        end
    end
    
end

