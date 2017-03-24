function [V,I,t,spike_times,firing_rate,spikes_per_cycle] = ...
    LIFNeuron(amp,Fc,dt,dur,Vrest,Vthresh,Vspike,Vreset,tau)

    % define input current as Fc Hz sinusoid with amplitude amp
    t = 0:dt:dur; % time vector
    I = amp*sin(2*pi*Fc*t); % injected current
    I(I < 0) = 0; % set negative values to 0

    % set the initial voltage
    V = zeros(1,length(t)); % voltage vector
    V(1) = Vrest;

    spike_count = 0; % spike count
    spike_times = []; 
    
    i = 2;
    while i <= length(t)
        % compute the change in voltage
        dV = (-(V(i-1) - Vrest) + I(i))*(dt/tau);
        V(i) = V(i-1) + dV;
        
        if V(i) > Vthresh % cell spiked
            spike_count = spike_count + 1;
            V(i) = Vspike;
            spike_times = [spike_times t(i)]; % record the spike time
            if i ~= length(t)
                V(i+1) = Vreset; % reset the voltage
                i = i + 1; % skip the next time step
            end
        end
        i = i + 1;
    end
    
    firing_rate = spike_count/dur; % firing rate
    spikes_per_cycle = spike_count/(Fc*dur); % spikes per cycle
     
    % convert to ms
    spike_times = spike_times*1000;
    t = t*1000;
    
end
