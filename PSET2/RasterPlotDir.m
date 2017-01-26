function RasterPlotDir(spikes,ndir,dir)

% plots a raster for a cell array of spike times

% INPUTS:
% spikes - cell array of spike times
% ndir - number of directions
% dir - column vector holding direction of each trial
    
    colors = hsv(ndir); % color matrix

    hold on; 
    
    num_trials = length(spikes);
    for i = 1:num_trials
        trial = spikes{i};
        len_trial = length(trial);
        color = colors((dir(i)/22.5)+1,:);
        
        for j = 1:len_trial
           spkx = [trial(j),trial(j)];
           spky = [i,i+0.9];
           line(spkx,spky,'color',color,'LineWidth',2);
        end
    end
   
    % plot horizontal line to separate directions
    direction_lines = 6:5:76;

    for i = 1:length(direction_lines)
        line([0,0.45],[direction_lines(i),direction_lines(i)],'color','k','LineWidth',0.25);
    end
    
    
end

