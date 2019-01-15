function numErrors(sample)

% Written by Dianna Amasino
% Finds the number of errors per subject, with errors defined as trials in 
% which choice went against subject value based on fitted k 

%Input: sample = 1 for primary sample, 2 for replication, behavioral data,
%subjective value data, discount rate output
%Output: file with the number of errors per subject

dataPath=pwd; %adapt to your location
cd(dataPath)
if sample ==1 %Primary sample
    load('amasinoETAl_behavior.csv') %load primary sample data
    data=amasinoETAl_behavior;
    load('subjVals.csv')
    load('allLogk.csv')
else % replication sample
    load('amasinoETAl_behavior_rep.csv')
    data=amasinoETAl_behavior_rep;
    load('subjVals_rep.csv')
    subjVals=subjVals_rep;
    load('allLogk_rep.csv')
    allLogk=allLogk_rep;
end

subj=1:data(end,1);
for i = 1:length(subj)  %loop over all subjects
    if ~isnan(allLogk(subj(i)))
        sub=find(data(:,1)==subj(i));       
        sv=subjVals(sub,2)>subjVals(sub,3); %If LL>SS = 1, if SS>LL = 0
        errors=find(sv~=data(sub,6) & ~isnan(data(sub,7)) & data(sub,2)<10); %exclude non-responses, trial where immed opt=10   
        %Find total # of errors
        numErrs(i)=length(errors);
        
    else
        numErrs(i)=NaN;
    end
end  

if sample==1
    csvwrite('numErrs.csv',numErrs')
else
    csvwrite('numErrs_rep.csv',numErrs')
end