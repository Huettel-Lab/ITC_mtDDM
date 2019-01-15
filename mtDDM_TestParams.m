function output = mtDDM_TestParams(data,params,RTBins,respOpts)
%% simulate ddm for each value difference

% Pre-allocate matrix sizes
FreqMatrix = zeros(length(respOpts),size(RTBins,2) - 1,length(params.possDiffs));
choice = zeros(params.nSims, length(params.possDiffs));
rt = zeros(params.nSims, length(params.possDiffs));

for d = 1:length(params.possDiffs) %Loop over all experienced values
   
    % Simulate choices and response times
    if params.type==0 %attribute-wise simulations
        [choice(:,d), rt(:,d)] = mtDDM_DriftSim(params,params.possDiffs(d,1),...
        params.possDiffs(d,2)); % Get simulated choices and RT for each value difference
    
    else %option-wise simulations
        [choice(:,d), rt(:,d)] = mtDDM_DriftSim_opt(params,params.possDiffs(d,1),...
            params.possDiffs(d,2), params.possDiffs(d,3),params.possDiffs(d,4),...
            data.avgA,data.avgT);
    end
    
    [~,whichRTBin]=histc(rt(:,d),RTBins(1,:)); % Find each simulated RT's RT bin
    whichChoiceBin = (choice(:,d)+3)/2; %Transform choices from [-1,1] to [1,2]   
    for r = 1:length(respOpts) % for each response option 
        for t = 1:max(whichRTBin)  %For each rt bin
        % count frequency of each response type within each RT bin from simulations
            FreqMatrix(r,t,d) = sum(whichChoiceBin == r & whichRTBin==t);
        end
    end  
end
% Normalize frequency bins to proportion of responses in each (RT X choice) bin
% by dividing by the number of experienced value differences x # simulations 
propMatrix = FreqMatrix./ (length(params.possDiffs)*params.nSims);
output.choice=choice;
output.rt=rt;
output.FreqMatrixNew=FreqMatrix;
output.propMatrix=propMatrix;

%% Calculate the log likelihood for each trial 

NLLDataProb = zeros(length(data.resp),1); % preallocate vector for NLL for each trial
alltheRespBins=(data.resp+3)/2; % convert response from [-1,1] to [1,2]
[~,allTheRTBins]=histc(data.rt,RTBins(1,:)); % Find each RT's RT bin (actual data)

for trial = 1:length(data.resp) % for each actual trial
    Respbin = alltheRespBins(trial); % Find that trial's reponse 
    RTbin = allTheRTBins(trial); % Find that trial's RT bin
    if params.type==0 %attribute-wise
        valBin = find(params.possDiffs(:,1) == data.DA(trial) & params.possDiffs(:,2) ...
        == data.DT(trial)); % Find that trial's value difference for amount and time 
    else % option-wise
        valBin = find(params.possDiffs(:,1) == data.AL(trial) & params.possDiffs(:,2) == data.AR(trial)...
        & params.possDiffs(:,3) == data.TL(trial) & params.possDiffs(:,4) == data.TR(trial));
    end
    
    % Assign the probability by finding the log of probability based on the 
    % choice/rt bin and value difference. The higher those are for each choice, 
    % the lower the NLL. Just sampling the experienced trials.
   if Respbin>0 && RTbin>0 && any(valBin)
       NLLDataProb(trial) = -log(propMatrix(Respbin,RTbin,valBin) + eps); % eps is to avoid infinite log value
   end
end

% sum all negative log likelihoods across trials for that set of parameters
output.negLL = sum(NLLDataProb);

end