function ET_errorAnalysis(sample)

% Written by Dianna Amasino
% Finds response time, and eye tracking index metrics for different types of 
% errors and corresponding, similar correct trials
% Error trials are trials where choice went against subjective value based on fitted k

% Inputs: sample = 1 for primary sample, 2 for replication sample,
% behavioral data, subjective value data, discount rate data, and eye
% tracking indices data
% Outputs: A matrix of average response times and eye tracking indices for
% error and correct trials. The order of columns is:
    % LL error response time, SS correct response time;
    % LL error Option Index, SS correct Option index;
    % LL error Attribute Index, SS correct Attribute index;
    % LL error Payne Index, SS correct Payne Index
    % SS error response time, LL correct response time;
    % SS error Option Index, LL correct Option index;
    % SS error Attribute Index, LL correct Attribute index;
    % SS error Payne Index, LL correct Payne Index

% Only include subjects with eye tracking data and fittable k-vals
dataPath=pwd; %adapt to your location
cd(dataPath)
if sample ==1 %Primary sample
    load('amasinoETAl_behavior.csv') %load primary sample data
    data=amasinoETAl_behavior;
    load('subjVals.csv')
    load('allLogk.csv')
    load('ET_indices.csv')
    subj=[1:5,7:15,17,19,21,22,24:44,47:49,51:56,58,60:63,68,70:73,75:77,80:85,87:91,93:98,100,103:109,111:117];
else % replication sample
    load('amasinoETAl_behavior_rep.csv')
    data=amasinoETAl_behavior_rep;
    load('subjVals_rep.csv')
    subjVals=subjVals_rep;
    load('allLogk_rep.csv')
    allLogk=allLogk_rep;
    load('ET_indices_rep.csv')
    ET_indices=ET_indices_rep;
    subj=[1,5:9,11:13,16,17,19,21:23,26:38,40,43,44,46:49,52,53,58,60,62:66,68:71,73,74,76:81,83:85,87:90,92,93,96,97,99];
end
       
for i = 1:length(subj)  %loop over all subjects
    sub=find(data(:,1)==subj(i));
    sub2=find(ET_indices(:,1)==subj(i)); %eye tracking subject values

    sv=subjVals(sub,2)>subjVals(sub,3); %If LL>SS = 1, if SS>LL = 0
    errors=find(sv~=data(sub,6) & ~isnan(data(sub,7)) & data(sub,2)<10); %exclude non-responses, trial where immed opt=10   
    correct=find(sv==data(sub,6) & ~isnan(data(sub,7)) & data(sub,2)<10);
    %Find error type and correct choice that is most similar in terms
    %of difference in sujbective value
    errorType{i}=sv(errors)-data(sub(errors),6); % +1 means impatient error, -1 means patient error
    errorSize{i}=subjVals(sub(errors),2)-subjVals(sub(errors),3); % find distance between correct resp and error 
    corrSize{i}=subjVals(sub(correct),2)-subjVals(sub(correct),3); % + means chose delay, - means chose immed

    corInd=[]; %gives non-error trials most similar to error trials by difference in SV
    for j=1:length(errorSize{i}) %Loop over each error
        %find first min dist in SV non-error
        if errorType{i}(j)>0 %impatient error, chose SS when LL>SS
            opp_resps=find(data(sub(correct),6)==1); %Chose LL correctly
        else %patient error
            opp_resps=find(data(sub(correct),6)==0); %Chose SS correctly
        end
        %Find correct trial with minimum difference in SV compared to error trial
        if ~isempty(opp_resps) %if there's at least 1 corresponding correct trial
            corInd(j)=find(abs(corrSize{i}(opp_resps)-errorSize{i}(j))==min(abs(corrSize{i}(opp_resps)-errorSize{i}(j))),1);
            corInd(j)=opp_resps(corInd(j));
        else
            corInd(j)=0;
        end
    end

    %Negative error means a patient error (chose LL when SS was correct), whereas
    %positive error means an impatient error (chose SS when LL was correct)
    negErrInd=find(errorType{i}<0);       
    if length(negErrInd)>2 % At least 3 of a given type of error
        negInd=sub(errors(negErrInd));
        negIndET=sub2(errors(negErrInd));
        negIndCor=sub(correct(corInd(negErrInd)));
        negIndCorET=sub2(correct(corInd(negErrInd)));
        %RT for error, RT for corresponding correct trial
        errorMatrix(i,1:2)=[mean(data(negInd,7)) mean(data(negIndCor,7))];
        if sum(~isnan(ET_indices(negIndET,3)))>2 %make sure ET data for at least 3 errors exists
            %Option Ind for error, correct; Attribute ind for error, correct
            errorMatrix(i,3:6)=[nanmean(ET_indices(negIndET,3)) nanmean(ET_indices(negIndCorET,3)) ...
                nanmean(ET_indices(negIndET,4)) nanmean(ET_indices(negIndCorET,4))];
        else %Not enough data for errors
            errorMatrix(i,3:6)=[NaN NaN NaN NaN];
        end

        if sum(~isnan(ET_indices(negIndET,5)))>2 %make sure PI data for at least 3 errors exists
            errorMatrix(i,7:8)=[nanmean(ET_indices(negIndET,5)) nanmean(ET_indices(negIndCorET,5))];
        else %Not enough data for errors
            errorMatrix(i,7:8)=[NaN NaN];
        end
    else %Not enough errors
        errorMatrix(i,1:8)=[NaN NaN NaN NaN NaN NaN NaN NaN];
    end

    posErrInd=find(errorType{i}>0);
    if length(posErrInd)>2 % At least 3 of a given type of error
        posInd=sub(errors(posErrInd));
        posIndET=sub2(errors(posErrInd));
        posIndCor=sub(correct(corInd(posErrInd)));
        posIndCorET=sub2(correct(corInd(posErrInd)));
        %RT for error, RT for corresponding correct trial
        errorMatrix(i,9:10)=[mean(data(posInd,7)) mean(data(posIndCor,7))];
        if sum(~isnan(ET_indices(posIndET,3)))>2 %make sure ET data for at least 3 errors exists
            %Option Ind for error, correct; Attribute ind for error, correct
            errorMatrix(i,11:14)=[nanmean(ET_indices(posIndET,3)) nanmean(ET_indices(posIndCorET,3)) ...
                nanmean(ET_indices(posIndET,4)) nanmean(ET_indices(posIndCorET,4))];
        else %Not enough data for errors
            errorMatrix(i,11:14)=[NaN NaN NaN NaN];
        end

        if sum(~isnan(ET_indices(posIndET,5)))>2 %make sure PI data for at least 3 errors exists
            errorMatrix(i,15:16)=[nanmean(ET_indices(posIndET,5)) nanmean(ET_indices(posIndCorET,5))];
        else %Not enough data for errors
            errorMatrix(i,15:16)=[NaN NaN];
        end
    else %Not enough errors
        errorMatrix(i,9:16)=[NaN NaN NaN NaN NaN NaN NaN NaN];
    end
end  

if sample==1
    csvwrite('errorMatrix.csv',errorMatrix)
else
    csvwrite('errorMatrix_rep.csv',errorMatrix)
end