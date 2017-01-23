function PSTH(spikes,binsize)

% plots a peri-stimulus time histogram based on the given bin size
% represents the time-varying "instantaneous firing rate"

% INPUTS:
% spikes - cell array of spike times
% binsize - bin size
% fs - fontsize

    % determine duration of observation - faster way?
    % note: some trials are empty
    duration = spikes{1}(end);
    ntrials = length(spikes);
    for i = 2:ntrials
        if ~isempty(spikes{i})   
            curr_duration = spikes{i}(end);
            if curr_duration > duration
                duration = curr_duration;
            end
        end
    end

    % create vector to hold spikes in each bin
    binspikes = zeros(1,ceil(duration/binsize));
    binranges = 0:binsize:duration;
    
    % place spikes of each trial in proper bin
    for i = 1:ntrials
        
        bincounts = histc(spikes{i},binranges);
        binspikes = binspikes + bincounts;

    end
    
    % create histogram with bar height of bin as numspikes/(numtrials*binsize)
    % to show instantaneous firing rate
    bar(binranges,binspikes/(ntrials*binsize),'histc');

    



