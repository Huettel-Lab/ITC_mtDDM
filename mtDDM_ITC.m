function mtDDM_ITC(subj,sample,type,resume)
% REVISION HISTORY
% Written by Nikki Sullivan
% Adapted by Dianna Amasino for intertemporal choice

% Fits a multi-attribute, multi-latency DDM:
% driftA - drift slope for amount information
% driftT - drift slope for time information
% latA - latency for amount information
% latT - latency for time information
% bounds - decision boundaries
% Bias is held constant at an assumed 0 (no bias for left/right fitting)
% Noise is Gaussian noise centered at 0 with a standard deviation held constant 
% at 0.1. Because the standard deviation is held constant, noise acts as a 
% scaling factor, and the drift slopes and bounds are fit relative to it.

% Uses a grid search procedure: 
% 1) Start with a coarse grid to get the parameter range for all subjects
% 2) Use a finer grid around the top fitting values for each subject

% Calls the mtDDM_TestParams and mtDDM_DriftSim or mtDDM_DrifftSim_opt scripts

% Inputs: 
% subj - subject ID to fit separately for each participant
% sample - 1 (primary) or 2 (replication) data
% type - 0 (attribute-wise) or 1 (option-wise) model form
% resume - 0 or empty if starting, 1 if resuming a previously stopped run
% The script loads the behavioral data and adapts it to the format needed
% for model fitting--it requires options seen as well as choice and
% response time.

% Outputs:
% Saves a file in a 'ddmOutputs' folder with the file name attDDM_[subj] or attDDMrep_[subj] 
% for attribute-wise DDMs for primary and replications sample, respectively.
% The file will be optDDM_[subj] or optDDMrep_[subj] for option-wise DDMs.
% If running multiple models on the same subject, the output file name should be changed.
% DDM parameters specified above
% Measures of fit:
% logL - Log likelihood
% AIC - Akaike Information Criterion
% BIC - Bayesian Information Criterion

%% Set seed for randomization
try
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
catch
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
end
%% Load and format subject specific data

dataPath=pwd; %adapt to your location
cd(dataPath)
if sample==1 %primary sample
    load('amasinoETAl_behavior.csv') %load primary sample data
    data=amasinoETAl_behavior;
else % sample==2
    load('amasinoETAl_behavior_rep.csv') %load replication data
    data=amasinoETAl_behavior_rep;
end

% Specify RT, choices, and options for a given subject
ind=data(:,1)==subj & ~isnan(data(:,6)); % Get the data for that subject excluding non-responses
ddmData.rt=data(ind,7); % Response times
side=data(ind,8); % side 1 = LL on left, 0 = LL on right; replication 1,3 = LL on top, 2,4 = LL on bottom
resp=data(ind,6); % resp 1 = LL, 0 = SS
if sample==1 %primary sample
    ddmData.resp(side==1 & resp==1)=1; % Chose left option (LL)
    ddmData.resp(side==1 & resp==0)=-1; % Chose right option (SS)
    ddmData.resp(side==0 & resp==1)=-1; % Chose right option (LL)
    ddmData.resp(side==0 & resp==0)=1; % Chose left option (SS)
    ddmData.AL=(side-1).*-data(ind,2); % When side = 0, SS on left
    ddmData.AL(ddmData.AL==0)=10; %Everything else LL on left
    ddmData.AR=side.*data(ind,2);% When side = 1, SS on right
    ddmData.AR(ddmData.AR==0)=10; %Everything else LL on right
    ddmData.TL=side.*data(ind,5); %When side = 1, LL on left, otherwise 0 time delay for SS
    ddmData.TR=(side-1).*-data(ind,5); %When side = 0, LL on right otherwise 0 time delay for SS
else % sample==2 replication sample
    ddmData.resp((side==1 | side==3) & resp==1)=1; % Chose top option (LL)
    ddmData.resp((side==1 | side==3) & resp==0)=-1; % Chose bottom option (SS)
    ddmData.resp((side==2 | side==4) & resp==1)=-1; % Chose bottom option (LL)
    ddmData.resp((side==2 | side==4) & resp==0)=1; % Chose top option (SS)
    ddmData.AL=round((mod(side,2)-1).*-data(ind,2)); % When side = 2,4, SS on top (left)
    ddmData.AL(ddmData.AL==0)=10; %Everything else LL on top (left)
    ddmData.AR=round(mod(side,2).*data(ind,2));% When side = 1,3 SS on bottom (right)
    ddmData.AR(ddmData.AR==0)=10; %Everything else LL on right
    ddmData.TL=mod(side,2).*data(ind,5); %When side = 1,3 LL on top (left), otherwise 0 time delay for SS
    ddmData.TR=(mod(side,2)-1).*-data(ind,5); %When side ind = 2,4 LL on bottom (right) otherwise 0 time delay for SS
end

% Find mean and RT quantiles
meanRT=mean(ddmData.rt);
lowQ=quantile(ddmData.rt,.025);
highQ=quantile(ddmData.rt,.975);
quants=100*quantile(ddmData.rt,linspace(.025,.975,8)); %Put quantiles in centiseconds
index=find(ddmData.rt>lowQ & ddmData.rt<highQ); %find 95% interval to exclude extreme RTs

%The amount differences are rounded below to reduce the # of simulations and simplify, 
%but it's not necessary and you might not want to do so with fewer combinations of options

%Format data to exclude extreme RTs and be in correct format
ddmData.rt=ddmData.rt(index)*100;% put in centiseconds since that's the step size
%Dif Amount (L-R) normalized to (-1,1), with 9.5 being maximum difference
ddmData.DA=round((ddmData.AL(index)-ddmData.AR(index))/9.5,1); %Round to just get $1 differences 
%Dif time (L-R) normalized to (-1,1), with 365 being maximum difference 
ddmData.DT=-(ddmData.TL(index)-ddmData.TR(index))/365; % negative because larger time is worse
ddmData.resp=ddmData.resp(index); %response input left = 1, right = -1
ddmData.AL=ddmData.AL(index); %Amount Left
ddmData.AR=ddmData.AR(index);%Amount right
ddmData.TL=ddmData.TL(index); %Time left
ddmData.TR=ddmData.TR(index); % Time right
ddmData.avgA=mean([ddmData.AL; ddmData.AR]);
ddmData.avgT=mean([ddmData.TL; ddmData.TR]);

%% Set DDM parameter ranges to test

% Check to see if restarting script to pick up where it left off
if exist('resume','var') 
    if resume
        if sample==1 && type==0;
            if exist(fullfile(dataPath,'ddmOutputs', ['attDDM_' num2str(subj) '.mat']),'file')
                load(fullfile(dataPath,'ddmOutputs',['attDDM_' num2str(subj)]))
            else
                resume=false;
            end
        elseif sample==1 && type==1;
            if exist(fullfile(dataPath,'ddmOutputs', ['optDDM_' num2str(subj) '.mat']),'file')
                load(fullfile(dataPath,'ddmOutputs',['optDDM_' num2str(subj)]))
            else
                resume=false;
            end
        elseif sample==2 && type==0;
            if exist(fullfile(dataPath,'ddmOutputs', ['attDDMrep_' num2str(subj) '.mat']),'file')
                load(fullfile(dataPath,'ddmOutputs',['attDDMrep_' num2str(subj)]))
            else
                resume=false;
            end
        elseif sample==2 && type==1;
            if exist(fullfile(dataPath,'ddmOutputs', ['optDDMrep_' num2str(subj) '.mat']),'file')
                load(fullfile(dataPath,'ddmOutputs',['optDDMrep_' num2str(subj)]))
            else
                resume=false;
            end     
        end
    else
        resume=false;
    end

params.nSims = 1000; % Set # of simulations per value difference 
if type==0 % attribute-wise
    params.possDiffs = unique([ddmData.DA'; ddmData.DT']','rows'); % Set value differences experienced by participant
    params.type=0;
else % option-wise
    params.possDiffs = unique([ddmData.AL'; ddmData.AR'; ddmData.TL'; ddmData.TR']','rows'); %Set value combos experienced
    params.type=1;
end

% ddm parameter ranges to search over, coarse search can be ~10 for each
if type == 0 % attribute-wise
    driftSlopeAmt = linspace(.0001,.1,10); % drift slope per 10 ms
    driftSlopeTime = linspace(.0001,.1,10); % drift slope per 10 ms
else % option-wise
    driftSlopeAmt = logspace(-4,-.8,10); % drift slope per 10 ms
    driftSlopeTime = logspace(-4,-.8,10); % drift slope per 10 ms
end

latencyAmt = round(linspace(10,100*meanRT,10)); % 10 ms steps, start at 100 ms, max at mean RT
latencyTime = round(linspace(10,100*meanRT,10)); % 10 ms steps, start at 100 ms, max at mean RT 
bounds = linspace(1,3,10); % boundaries (assume symmetrical)

params.bias = 0; % starting bias.
params.rvsNoise = .1; % RVS noise standard deviation set to .1 per 10 ms time step 

% response & response time bins
RTBins=quants; %Use quantiles from above to define RT bins
respOpts = [-1, 1]; % [right, left] 
%% simulate each ddm & find best fitting one.

% Set all possible ddm parameter combinations in a matrix
allCombos=combvec(driftSlopeAmt,driftSlopeTime,latencyAmt,latencyTime,bounds)';

% Preallocate response variables
if ~resume
    negLL = NaN(length(allCombos),1); %NaN matrix to fill in
    out = cell(length(allCombos),1); %Output 
end

% Find the indices that haven't been run yet:
pairsNotRun = find(isnan(negLL)); % number of runs left

% Set frequency of saving (takes time so don't do it each run)
percDone=100*((1:length(pairsNotRun))/length(pairsNotRun));
xPerc=find(percDone>=1,1,'first'); % Save every 1 percent (can change to percDone>=5 or other value if desired).
saveFile=false(length(pairsNotRun),1); %don't save if not x percent
saveFile(1:xPerc:end)=true; %save every x percent
saveFile(end)=true; %save at end

for indRun = 1:length(pairsNotRun) %loop over everything that hasn't been run
    nC = pairsNotRun(indRun); %nC is the row index of parameter values to use in a given simulation
    
    params.d1 = allCombos(nC,1); % Amt drift slope
    params.d2 = allCombos(nC,2); % Time drift slope
    params.ndt1 = allCombos(nC,3); % Amt latency
    params.ndt2 = allCombos(nC,4); % Time latency
    params.bound = allCombos(nC,5); % Bounds (symmetrical)
   
    % now run it
    out{nC} = mtDDM_TestParams(ddmData,params,RTBins,respOpts); %Run mtDDM_TestParams script
    negLL(nC) = out{nC}.negLL; %assign values to negLL once it has been run
    
    if saveFile(indRun) %save every x percent as specified above
       mkdir(fullfile(dataPath,'ddmOutputs'))
        if sample==1;
            if type==0;
                save(fullfile(dataPath,'ddmOutputs',['attDDM_' num2str(subj)]),'negLL','allCombos');
            else
                save(fullfile(dataPath,'ddmOutputs',['optDDM_' num2str(subj)]),'negLL','allCombos');
            end
        else
            if type==0;
                save(fullfile(dataPath,'ddmOutputs',['attDDMrep_' num2str(subj)]),'negLL','allCombos');
            else
                save(fullfile(dataPath,'ddmOutputs',['optDDMrep_' num2str(subj)]),'negLL','allCombos');
            end
        end
    end
end
%% find best simulation

%Find minimum negative log likelihood
bestInd = find(negLL == nanmin(negLL));

%Find best parameter combinations
ddmOut.driftA = allCombos(bestInd,1);
ddmOut.driftT = allCombos(bestInd,2);
ddmOut.latA = allCombos(bestInd,3);
ddmOut.latT = allCombos(bestInd,4);
ddmOut.bounds = allCombos(bestInd,5);

%Calculate fit indices
nParams=5; % 2 drift slopes, 2 latencies, 1 bound (symmetrical)
ddmOut.logL = -negLL(bestInd);
ddmOut.G = 2*negLL(bestInd);
ddmOut.BIC = ddmOut.G + nParams*log(length(ddmData.resp));
ddmOut.AIC = ddmOut.G + 2*nParams;

% save the fitting procedure
ddmOut.ranges.bias = params.bias;
ddmOut.ranges.driftSlopeA = driftSlopeAmt;
ddmOut.ranges.driftSlopeT = driftSlopeTime;
ddmOut.ranges.latencyA = latencyAmt;
ddmOut.ranges.latencyT = latencyTime;
ddmOut.ranges.bounds = bounds;
ddmOut.finished = true;

%Save everything
 if sample==1;
    if type==0;
        save(fullfile(dataPath,'ddmOutputs',['attDDM_' num2str(subj)]),'ddmOut','negLL','allCombos');
    else
        save(fullfile(dataPath,'ddmOutputs',['optDDM_' num2str(subj)]),'ddmOut','negLL','allCombos');
    end
else
    if type==0;
        save(fullfile(dataPath,'ddmOutputs',['attDDMrep_' num2str(subj)]),'ddmOut','negLL','allCombos');
    else
        save(fullfile(dataPath,'ddmOutputs',['optDDMrep_' num2str(subj)]),'ddmOut','negLL','allCombos');
    end
 end

end