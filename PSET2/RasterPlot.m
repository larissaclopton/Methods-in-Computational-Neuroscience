function RasterPlot(spikes,color)

% plots a raster for a cell array of spike times

% INPUTS:
% spikes - cell array of spike times
% color - the color for the plot
   
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
    
    % draw lines to separate pre-trail, sample, delay and test periods
    line([0 0],[1 num_trials+1],'color','k','LineWidth',3);
    line([0.65 0.65],[1 num_trials+1],'color','k','LineWidth',3);
    line([1.65 1.65],[1 num_trials+1],'color','k','LineWidth',3);
    axis([-0.5 2 1 num_trials]);
end
