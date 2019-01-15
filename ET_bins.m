function ET_bins(sample)
%Written by Dianna Amasino
% Splits up eye tracking data into time bins and finds the proportion of time
%looking left (or top) within those time bins split by the option chosen. 
%Also finds proportion of last fixations to the left (or top) given the side
%that was chosen.

%Inputs: sample (1 is primary sample, 2 is replication sample), behavioral data,
%and corrected eye tracking data
%Outputs: Eye tracking bins with average proportion of looking left for
%each subject over trials when they chose the left option (binsL.csv) or
%the right option (binsR.csv). Also, the proportion of left choice based on
%the last fixation being to the left (looksL.csv) or right (looksR.csv)

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

%AOIs 1, 2 are on the left ; 3, 4 are the right 
%AOIs 1,3 are on the top, 2 & 4 are on the bottom
if sample==1
    gazeSide(data(:,6)==1 | data(:,6)==2)=1; %code left AOI gazes as 1, right gazes as 0
    gazeSide(data(:,6)==3 | data(:,6)==4)=0; %code left AOI gazes as 1, right gazes as 0
elseif sample==2
    gazeSide(data(:,6)==1 | data(:,6)==3)=1; %code left AOI gazes as 1, right gazes as 0
    gazeSide(data(:,6)==2 | data(:,6)==4)=0; %code left AOI gazes as 1, right gazes as 0
end

%Loop over each subject
for i = 1:length(subj)
    choseLeft=[];
    choseRight=[];
    lastLookL=[];
    lastLookR=[];
    sub = find(data(:,1)==subj(i)); %Find subject-specific data points
    %Loop over each trial
    for j=1:trials;
        row=j+(subj(i)-1)*141; %account for continuous rows in beh data
        %row2=j+(i-1)*141;
        ind=find(data(sub,2)==j);
        
        %Split trial into fifths, look at how gaze evolves over time
        if length(ind)>4 %at least 5 gaze points
            fifth=round(length(ind)/5); 
            fifths=[nanmean(gazeSide(sub(ind(1:fifth)))) nanmean(gazeSide(sub(ind(fifth+1:2*fifth))))...
                nanmean(gazeSide(sub(ind(2*fifth+1:3*fifth)))) nanmean(gazeSide(sub(ind(3*fifth+1:4*fifth))))...
                nanmean(gazeSide(sub(ind(4*fifth+1:end))))];
            
            if beh(row,9)==1 %chose left  
                choseLeft=[choseLeft; fifths]; % average left looking time for each time point
            elseif beh(row,9)==0 %chose right
                choseRight=[choseRight; fifths];
            end
            
            % Look at % left choices based on direction of last fixation
            if isempty(gazeSide(sub(ind))) || isnan(gazeSide(sub(ind(end)))) %if look is empty or the end is NaN
                % do nothing
            elseif gazeSide(sub(ind(end)))==0 %last fixation was right
                lastLookR=[lastLookR; beh(row,9)];
            elseif gazeSide(sub(ind(end)))==1 %last fixation was left
                lastLookL=[lastLookL; beh(row,9)];
            end
        end
    end
    binsL(i,:)=[subj(i) nanmean(choseLeft)]; 
    binsR(i,:)=[subj(i) nanmean(choseRight)];
    looksL(i,:)=[subj(i) mean(lastLookL)];
    looksR(i,:)=[subj(i) mean(lastLookR)];
end
if sample==1;
    csvwrite('binsL.csv',binsL) 
    csvwrite('binsR.csv',binsR)
    csvwrite('looksL.csv',looksL)
    csvwrite('looksR.csv',looksR)
else
    csvwrite('binsL_rep.csv',binsL) 
    csvwrite('binsR_rep.csv',binsR)
    csvwrite('looksL_rep.csv',looksL)
    csvwrite('looksR_rep.csv',looksR)
end