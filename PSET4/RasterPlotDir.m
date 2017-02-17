function RasterPlotDir(spikes,times,ndir,ntrials)

% INPUTS:
% spikes - 3D matrix of boolean spikes (time x ndir x trials)
% times - time of each bin
% ndir - number of directions
% ntrials - number of trials
    
    colors = hsv(ndir); % color matrix
    len_trial = length(times);

    hold on; 
     
    % sort by stimulus direction
    for i = 1:ndir
        color = colors(i,:);
        for j = 1:ntrials
            for k = 1:len_trial
                if spikes(k,i,j)
                    spkx = [1000*times(k),1000*times(k)];
                    spky = [(i-1)*ntrials+(j-1),(i-1)*ntrials+(j-1)+0.9];
                    line(spkx,spky,'color',color,'LineWidth',5);
                end
            end
        end
    end
   
    % plot horizontal line to separate directions
    direction_lines = 184:184:13*184;

    for i = 1:length(direction_lines)
        line([0,1000*times(end)],[direction_lines(i),direction_lines(i)],'color','k','LineWidth',0.25);
    end
    
    axis([0,1000*times(end),0,13*184]);
    
end

