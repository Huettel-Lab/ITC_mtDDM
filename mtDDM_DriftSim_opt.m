function [choices rts] = mtDDM_DriftSim_opt(params,A1,A2,T1,T2,avgA,avgT)
% Run ddm simulations using given drift parameters
% Equation of interest assumes option-wise comparison, hyperbolic combination of
% amount and time: RVS(t) = RVS(t-1) + dA*A(l)/[1+dT*T(l)] + dA*A(r)/[1+dT*T(r)] 
% where actual A and T don't come online until a latency is passed for each, 
% separately; A and T have opposite values, and noise is added at every time point
% Note for option-wise discounting, you can't have amount or time as 0
% because they are dependent on each other in the equation, so we use the
% average amount and time (across the experiment) for pre-latency accumulation.
%% initialize/preallocate

nSteps = 1001; % time steps--10 ms, account for up to 10 s

%Take a random number for each time step for each simulation and multiply
%by noise standard deviation (0.1). To use in simulations below
draws = randn(params.nSims,nSteps).*(params.rvsNoise); % noise

%% run all trials simultaneously

% create path without noise
slope=[];
if params.ndt1 < params.ndt2 %latency amt < latency time
    % define path as zeros before amount latency (no info accumulation),
    % amount drift rate and value plus time drift rate and average time value 
    % between amount latency and time latency, and both drift rates and actual 
    % trial values for the rest of the trial. 
    slope = [slope zeros(1, params.ndt1) repmat((params.d1*A1/(1+params.d2*avgT))-...
        (params.d1*A2/(1+params.d2*avgT)), 1, params.ndt2-params.ndt1)]; 
    slope = [slope repmat((params.d1*A1/(1+params.d2*T1))-...
        (params.d1*A2/(1+params.d2*T2)), 1, nSteps - params.ndt2)];
elseif params.ndt1 > params.ndt2 % Latency amt > latency time
    % define path as zeros before time latency (no info accumulation),
    % time drift rate and value plus amount drift rate and average amount value 
    % between time latency and amount latency, and both drift rates and actual 
    % trial values for the rest of the trial. 
    slope = [slope zeros(1,params.ndt2) repmat((params.d1*avgA/(1+params.d2*T1))-...
        (params.d1*avgA/(1+params.d2*T2)),1, params.ndt1 - params.ndt2)];
    slope = [slope repmat((params.d1*A1/(1+params.d2*T1))-...
        (params.d1*A2/(1+params.d2*T2)), 1, nSteps - params.ndt1)];
else % Latencies are equal so define path as zeros before drift rates and both drift rates after latency
    slope = [slope zeros(1,params.ndt2) repmat((params.d1*A1/(1+params.d2*T1))-...
        (params.d1*A2/(1+params.d2*T2)),1, nSteps - params.ndt2)];
end

% add random noise to each path and get cumulative time points for each path
rvs = cumsum(draws + repmat(slope,params.nSims,1),2);

% get time at which each path crosses params.barrier
[didCross, ind] = max(abs(rvs)>=params.bound,[],2); % did relative value signal cross boundary/threshold and when?
rts = didCross .* ind; % remove paths that didn't cross params.barrier, ind is when it crossed

% find choice made by each path
choices = zeros(size(rts)); %making vector
pathSigns = sign(rvs); %is each path +/-?
for i = find(rts>0)' %only rts that finished
    choices(i) = pathSigns(i,rts(i)); %assign value to choice [-1,1]
end

end