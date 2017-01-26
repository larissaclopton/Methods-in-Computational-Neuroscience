function NeurometricFunc(pref,stimuli,spikes)

% computes and plots the neurometric function

% INPUTS:
% pref - the preferred direction or orientation of the neuron
% stimuli - a column vector of stimulus conditions
% spikes - a cell array of spike times

    % compute angular differences with respect to preferred
    ang_diffs = abs(stimuli - pref);
    ang_diff = unique(ang_diffs); % the unique elements of ang_diffs
    ndiff = length(ang_diff); % number of angular differences

    pref_samples = spikes(ang_diffs == 0);
    npref = length(pref_samples); % number of preferred trials
    p_correct = zeros(1,ndiff-1); % proportion correct for each angular difference

    % isolate sample for each angular difference
    for i = 2:ndiff
        nonpref_samples = spikes(ang_diffs == ang_diff(i));
        nnonpref = length(nonpref_samples); % number of nonpreferred trials
        correct = 0;
        % compare all combinations between pref_samples and non_pref samples
        % correct if more spikes in pref_sample than nonpref_sample
        for j = 1:npref
            for k = 1:nnonpref
                if length(pref_samples{j}) > length(nonpref_samples{k})
                    correct = correct + 1;
                end
            end
        end
        p_correct(i-1) = correct/(npref*nnonpref);
    end

    % plot this as proportion correct against angular difference
    plot(ang_diff(2:end),p_correct);

end