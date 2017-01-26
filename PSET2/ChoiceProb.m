function ChoiceProb(pref,nonpref)

% returns the choice probability for each period in a trial

% INPUTS:
% pref - cell array of spike times for preferred choice
% nonpref - cell array of spike times for nonpreferred choice

    period = {['pre-trial'] ['sample'] ['delay'] ['test']};
    window = {[-0.5 0] [0 0.65] [0.65 1.65] [1.65 2]}; % time windows for each period

    npref = length(pref);
    nnonpref = length(nonpref);

    for i = 1:4
        correct = 0;
        % compare all combinations of pref and nonpref pairs
        for j = 1:npref
            pref_spikes = find(pref{j} > window{i}(1) &...
                        pref{j} <= window{i}(2));
            for k = 1:nnonpref  
                nonpref_spikes = find(nonpref{k} > window{i}(1) &...
                        nonpref{k} <= window{i}(2));
                if length(pref_spikes) > length(nonpref_spikes)
                    correct = correct + 1;
                end  
            end
        end
        choice_prob = correct/(npref*nnonpref);
        fprintf(['Choice probability for ' period{i} ' period = ' num2str(choice_prob) '\n']);
    end

end