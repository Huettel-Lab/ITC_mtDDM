function subjectiveValue(sample)

% Written by Dianna Amasino
% Finds the trial by trial subject value of the smaller, sooner (SS) and 
% larger, later (LL) options and of the left and right options based on 
% each individual subject's fitted discount rate, k.
% Also splits trials into 10 bins by difference in subjective value (SV) between 
% the left and right options and finds the proportion of left responses,
% the response times, and the number of fixations for each SV bin.

%Input: sample =1 for primary sample, 2 for replication, behavioral data,
%discount rates, and number of fixations per trial per subject. 
%Output: A matrix (subjVals.csv) with the subject ID, SV of the LL and SS options and SV
%of the left and right options. Also, the proportion of left choices (SVpropL.csv),
%the average response time (SVrt.csv) and average number of fixations(SVfixNums.csv) 
%for SV left - SV right option value bins from [-10, 10].

dataPath=pwd; %adapt to your location
cd(dataPath)
    if sample ==1 %Primary sample
        load('amasinoETAl_behavior.csv') %load primary sample data
        load('allLogk.csv')
        load('numFixes.csv')
        data=amasinoETAl_behavior;
        logk=allLogk;
    else % replication sample
        load('amasinoETAl_behavior_rep.csv')
        load('allLogk_rep.csv')
        data=amasinoETAl_behavior_rep;
        logk=allLogk_rep;
        load('numFixes_rep.csv')
        numFixes=numFixes_rep;
    end
    
    subj=1:data(end,1);
for i = 1:length(subj)  %loop over all subjects
    k=exp(logk(i));
    rows=find(data(:,1)==subj(i));
    svLL(rows,1)=data(rows,3)./(1+k.*data(rows,5));
    svSS(rows,1)=data(rows,2)./(1+k.*data(rows,4));
    if sample==1; %Left of right side
        rows0=find(data(rows,8)==0); % side = 0 means SS on left, LL on right
        rows1=find(data(rows,8)==1); % side = 1 means LL on left, SS on right
    else %Top or bottom
        rows0=find(data(rows,8)==2 | data(rows,8)==4); % (side == 2 or 4) means SS on top, LL on bottom
        rows1=find(data(rows,8)==1 | data(rows,8)==3); % (side == 1 or 3) means LL on top, SS on bottom
    end
    svLT(rows(rows0),1)=svSS(rows(rows0)); % Subjective value left or Top
    svRB(rows(rows0),1)=svLL(rows(rows0)); % Subjective value right or bottom 
    svLT(rows(rows1),1)=svLL(rows(rows1)); % Subjective value left or Top
    svRB(rows(rows1),1)=svSS(rows(rows1)); % Subjective value right or bottom 
    dSV(rows,1)=svLT(rows)-svRB(rows); %Difference in subjective value (Left - right or Top - bottom)
    
    %Put values into 10 bins of SV left - SV right  bins 
    for j=1:10; %for range of SV left - SV right is (-10,10) in increments of 2   
        ind=find(dSV(rows)>=(-12+2*j) & dSV(rows)<(-10+2*j)); 
        propL(j)=nanmean(data(rows(ind),9)); %Find proportion of left choices
        rt(j)=nanmean(data(rows(ind),7)); %find RT
        if sum(numFixes(:,1)==subj(i)) %This person has eye tracking
            rowsET=find(numFixes(:,1)==subj(i)); 
            fixNums(j)=nanmean(numFixes(rowsET(ind),2)); %Find avg. # fixations 
        else
            fixNums(j)=NaN;
        end
    end
    SVpropL(i,:)=[i propL];
    SVrt(i,:)=[i rt];
    SVfixNums(i,:)=[i fixNums];
end    
subjVals=[data(:,1) svLL svSS svLT svRB dSV];

 if sample==1;
    csvwrite('subjVals.csv',subjVals) % write log subjective values to csv
    csvwrite('SVpropL.csv',SVpropL)
    csvwrite('SVrt.csv',SVrt)
    csvwrite('SVfixNums.csv',SVfixNums)
else
    csvwrite('subjVals_rep.csv',subjVals) % write log subjective values to csv
    csvwrite('SVpropL_rep.csv',SVpropL)
    csvwrite('SVrt_rep.csv',SVrt)
    csvwrite('SVfixNums_rep.csv',SVfixNums)
end