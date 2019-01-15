function [choices rts] = mtDDM_DriftSim(params,v1,v2)
% Run ddm simulations using given drift parameters
% Equation of interest assumes attribute-wise comparison, linear combination 
% RVS(t)= RVS(t-1)+dA[A(l)-A(r)]+dT[T(l)-T(r)] where A and T don't come
% online until a latency is passed for each, separately, A and T have opposite values, 
% and noise is added at every time point

%% initialize/preallocate

nSteps = 1001; % time steps--10 ms, 1000 steps = 10s which is the maximum trial length
  
% Take a random number for each time step for each simulation and multiply
% by noise standard deviation (0.1). To use in simulations below
draws = randn(params.nSims,nSteps).*(params.rvsNoise); % noise

%% run all trials simultaneously

% create path without noise
slope=[];
if params.ndt1 < params.ndt2 %latency amt < latency time
    % define path as zeros before amount latency (no info accumulation),
    % amount drift rate between amount latency and time latency, and both 
    % drift rates for the rest of the trial
    slope = [slope zeros(1, params.ndt1) repmat(params.d1*v1, 1, params.ndt2-params.ndt1)]; 
    slope = [slope repmat((params.d1*v1+params.d2*v2), 1, nSteps - params.ndt2)];
elseif params.ndt1 > params.ndt2 % latency amt > latency time
    % define path as zeros before time latency (no info accumulation),
    % drift rate time between time latency and amt latency, and both drift 
    % rates for the rest of the trial
    slope = [slope zeros(1,params.ndt2) repmat(params.d2*v2,1, params.ndt1 - params.ndt2)];
    slope = [slope repmat((params.d1*v1+params.d2*v2), 1, nSteps - params.ndt1)];
else % Latencies are equal so define path as zeros before latency and both drift rates after latency
    slope = [slope zeros(1,params.ndt2) repmat((params.d1*v1+params.d2*v2),1, nSteps - params.ndt2)];
end

% add random noise to each path and get cumulative time points for each path
rvs = cumsum(draws + repmat(slope,params.nSims,1),2);

% get time at which each path crosses params.barrier
[didCross, ind] = max(abs(rvs)>=params.bound,[],2); % did relative value signal cross boundary/threshold and when?
rts = didCross .* ind; % remove paths that didn't cross params.barrier, ind is when it crossed

% find choice made by each path
choices = zeros(size(rts)); %making vector
pathSigns = sign(rvs); %is each path +/-?
for i = find(rts>0)'; %only rts that finished
    choices(i) = pathSigns(i,rts(i)); %assign value to choice [-1,1]
end

end