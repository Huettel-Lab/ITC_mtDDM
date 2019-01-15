function concatDDMs(sample,type)

% Written by Dianna Amasino
% Loops over each subject's DDM output and puts all outputs in one file

% Inputs: sample (1 for primary sample, 2 for replication sample),
% type of DDM (0 for attribute-wise, 1 for option-wise), DDM output
% Outputs: a file with all the subjects' DDM results

dataPath=pwd; %adapt to your location
cd(dataPath)
if sample ==1 %Primary sample
    subj=1:117;
else % replication sample
    subj=1:100;
end

for i = 1:length(subj)  %loop over all subjects 
    if sample==1;
        if type==0;
            load(fullfile(dataPath,'ddmOutputs',['attDDM_' num2str(subj(i)) '.mat']));
        else
            load(fullfile(dataPath,'ddmOutputs',['optDDM_' num2str(subj(i)) '.mat']));
        end
    else
        if type==0;
            load(fullfile(dataPath,'ddmOutputs',['attDDMrep_' num2str(subj(i)) '.mat']));
        else
            load(fullfile(dataPath,'ddmOutputs',['optDDMrep_' num2str(subj(i)) '.mat']));
        end
    end
    ddmOut(i,1:8)=[ddmOut.driftA ddmOut.driftT ddmOut.latA ddmOut.latT ddmOut.bounds ddmOut.logL ddmOut.BIC ddmOut.AIC];
end

%save data
if sample==1; % Primary sample
    if type==0; %attribute wise
        csvwrite('attDDM.csv',ddmOut)
    else %option-wise
        csvwrite('optDDM.csv',ddmOut)
    end
else % Replication sample
    if type==0; %attribute wise
        csvwrite('attDDM_rep.csv',ddmOut)
    else %option-wise
        csvwrite('optDDM_rep.csv',ddmOut)
    end
end