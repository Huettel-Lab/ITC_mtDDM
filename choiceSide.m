function choiceSide(sample)
% Written by Dianna Amasino
% Finds the average proportion of left choices for each subject

% Inputs: sample = 1 (primary) or 2 (replication), and behavioral data
% Outputs: a file with the average proportion of left choices for each subject

dataPath=pwd; %adapt to your location
    cd(dataPath)
    if sample ==1 %Primary sample
        load('amasinoETAl_behavior.csv') %load primary sample data
        data=amasinoETAl_behavior;
    else % replication sample
        load('amasinoETAl_behavior_rep.csv')
        data=amasinoETAl_behavior_rep;
    end
    
    subj=1:data(end,1);
    for i = 1:length(subj)  %loop over all subjects
        %find subject-specific data and exclude non-responses
        ind=data(:,1)==subj(i) & ~isnan(data(:,6));
        percL(i,1)=mean(data(ind,9));
    end
    
    if sample==1;
        csvwrite('percLeft.csv',percL) % write RTs to csv
    else
        csvwrite('percLeft_rep.csv',percL) % write RTs to csv
    end