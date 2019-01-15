function ET_indices(sample)
% Original version written by Dianna Amasino 2016, updated 2019

% This function takes in corrected/cleaned eye tracking data and calculates
% the Payne Index, Attribute Index, and Option Index.

% INPUT
% sample - 1 (primary) or 2 (replication) data
% col 1: subject number (matching behavioral data)
% col 2: Trial number
% col 3: x-position (pixels on screen)
% col 4: y-position (pixels on screen)
% col 5: fixation counter
% col 6: AOI number (1 for TL, 2 for BL, 3 for TR, 4 for BR)

% OUTPUT
% Trial by trial indices--Option index, Attribute index, and Payne index
% Average index per person
% Average % first fixations to the amount compared to time
% Average number of fixations per trial

%Screen variables
% screen dimension
% width_x = 1280;
% height_y = 1024; 

%Load data
dataPath=pwd; %adapt to your location
cd(dataPath)
if sample ==1 % Primary sample
    load('ET_adj.csv') 
    data=ET_adj;
    load('amasinoETAl_behavior.csv')
    beh=amasinoETAl_behavior;
    subj=[1:5,7:17,19,21:56,58,60:63,65,68,70:73,75:91,93:101,103:117]; 
else % Replication sample
    load('ET_adj_rep.csv')
    data=ET_adj_rep;
    load('amasinoETAl_behavior_rep.csv')
    beh=amasinoETAl_behavior_rep;
    subj=[1:3,5:9,11:17,19,21:24,26:38,40,41,43:49,51:55,58,60,62:71,73:81,83:94,96:100];
end

trials=141;
ind1=find(data(:,6)==1);
ind2=find(data(:,6)==2);
ind3=find(data(:,6)==3);
ind4=find(data(:,6)==4);
%loop over each subject
for i = 1:length(subj)
    sub = find(data(:,1)==subj(i)); %Find subject-specific data points
    %LL and SS switch sides trial to trial, so account for this in AOIs 
    % Make LLamt=1, LLtime=2, SSamt=4, SStime=8 in a 7th column
    for j=1:trials;
        row=j+(subj(i)-1)*141; %account for continuous rows in beh data
        row2=j+(i-1)*141;
        ind=find(data(sub,2)==j);
        
        if sample==1 %Primary sample has amounts on top, times on bottom
            if beh(row,8)==0; % side = 0 means SS on left, LL on right
                data(intersect(sub(ind),ind1),7)=4;
                data(intersect(sub(ind),ind2),7)=8;
                data(intersect(sub(ind),ind3),7)=1;
                data(intersect(sub(ind),ind4),7)=2;
            else % side = 1 means LL on left, SS on right
                data(intersect(sub(ind),ind1),7)=1;
                data(intersect(sub(ind),ind2),7)=2;
                data(intersect(sub(ind),ind3),7)=4;
                data(intersect(sub(ind),ind4),7)=8;
            end
        else %Replication had 4 options for sides
            if beh(row,8)==1; % side = 1 means LL on top, SS on bottom; amount on left and time on right
                data(intersect(sub(ind),ind1),7)=1;
                data(intersect(sub(ind),ind2),7)=4;
                data(intersect(sub(ind),ind3),7)=2;
                data(intersect(sub(ind),ind4),7)=8;
            elseif beh(row,8)==2; % side = 2 means SS on top, LL on bottom; amount on left and time on right
                data(intersect(sub(ind),ind1),7)=4;
                data(intersect(sub(ind),ind2),7)=1;
                data(intersect(sub(ind),ind3),7)=8;
                data(intersect(sub(ind),ind4),7)=2;
            elseif beh(row,8)==3; % side = 3 means LL on top, SS on bottom; time on left and amount on right
                data(intersect(sub(ind),ind1),7)=2;
                data(intersect(sub(ind),ind2),7)=8;
                data(intersect(sub(ind),ind3),7)=1;
                data(intersect(sub(ind),ind4),7)=4;
            else % side = 4 means SS on top, LL on bottom; time on left and amount on right
                data(intersect(sub(ind),ind1),7)=8;
                data(intersect(sub(ind),ind2),7)=2;
                data(intersect(sub(ind),ind3),7)=4;
                data(intersect(sub(ind),ind4),7)=1;
            end
        end
        data(sub(ind),8)=[0; diff(data(sub(ind),7))]; %Find transitions within each trial
        %Gaze measures
        gazeAmt=sum(data(sub(ind),7)==1)+sum(data(sub(ind),7)==4);
        gazeTime=sum(data(sub(ind),7)==2)+sum(data(sub(ind),7)==8);
        gazeLL=sum(data(sub(ind),7)==1)+sum(data(sub(ind),7)==2);
        gazeSS=sum(data(sub(ind),7)==4)+sum(data(sub(ind),7)==8);
        transitAtt=sum(abs(data(sub(ind),8))==3)+sum(abs(data(sub(ind),8))==6);
        transitOpt=sum(abs(data(sub(ind),8))==1)+sum(abs(data(sub(ind),8))==4);
        
        %Find subject-specific index for trial
        if ~isempty(ind)
            firstFix(row2,1:2)=[subj(i) data(sub(ind(1)),7)];
            numFix(row2,1:2)=[subj(i) 1+data(sub(ind(end)),5)-data(sub(ind(1)),5)];
        else
            firstFix(row2,1:2)=[subj(i) NaN];
            numFix(row2,1:2)=[subj(i) NaN];
        end
        OI_trial(row2,1)=(gazeSS-gazeLL)/(gazeSS+gazeLL);
        AI_trial(row2,1)=(gazeAmt-gazeTime)/(gazeAmt+gazeTime);
        PI_trial(row2,1)=(transitOpt-transitAtt)/(transitOpt+transitAtt);
        indices(row2,:)=[subj(i) row2-(i-1)*141 OI_trial(row2) AI_trial(row2) PI_trial(row2)];
    end
    % Average indices across trials 
    amountFixes(i)=sum(firstFix(indices(:,1)==subj(i),2)==1)+sum(firstFix(indices(:,1)==subj(i),2)==4);
    timeFixes(i)=sum(firstFix(indices(:,1)==subj(i),2)==2)+sum(firstFix(indices(:,1)==subj(i),2)==8);
    firstFixes(i,1)=amountFixes(i)/(amountFixes(i)+timeFixes(i));
    OI(i,1)= nanmean(indices(indices(:,1)==subj(i),3));
    AI(i,1)= nanmean(indices(indices(:,1)==subj(i),4));
    PI(i,1)= nanmean(indices(indices(:,1)==subj(i),5));        
end

%Save output
if sample==1
    csvwrite('ET_indices.csv',indices);
    csvwrite('OptInd.csv',OI);
    csvwrite('AttInd.csv',AI);
    csvwrite('PayneInd.csv',PI);
    csvwrite('firstFixes.csv',firstFixes)
    csvwrite('numFixes.csv',numFix)
else
    csvwrite('ET_indices_rep.csv',indices)
    csvwrite('OptInd_rep.csv',OI)
    csvwrite('AttInd_rep.csv',AI)
    csvwrite('PayneInd_rep.csv',PI)
    csvwrite('firstFixes_rep.csv',firstFixes)
    csvwrite('numFixes_rep.csv',numFix)
end
