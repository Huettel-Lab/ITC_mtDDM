function DDM_simulations(sample,type)

%Written by Dianna Amasino
% Runs simulations for each best-fitting DDM parameters for each subject and 
% finds the average proportion of LL choices as well as the response time
% quantiles. 

% Inputs: sample = 1 for primary, 2 for replication, type = 0 for
% attribute-wise DDM, 1 for option-wise DDM, behavioral data, DDM output
% Outputs: actual and simulated proportion of LL choices; actual and
% simulated RT quantiles for each subject

dataPath=pwd; %adapt to your location
cd(dataPath)

%Load data
if sample ==1 %Primary sample
    load('amasinoETAl_behavior.csv') %load primary sample data
    data=amasinoETAl_behavior;
    if type==0 %attribute-wise
        load('attDDM.csv')
        ddm=attDDM;
    else %option-wise
        load('optDDM.csv')
        ddm=optDDM;
    end
else % replication sample
    load('amasinoETAl_behavior_rep.csv')
    data=amasinoETAl_behavior_rep;
    if type==0 %attribute-wise
        load('attDDM_rep.csv')
        ddm=attDDM_rep;
    else %option-wise
        load('optDDM_rep.csv')
        ddm=optDDM_rep;
    end       
end
subj=1:data(end,1);

%Loop over each subject
for i = 1:length(subj)
    %Prep the data as in mtDDM_ITC.m
    % Specify RT, choices, and options for a given subject
    ind=data(:,1)==subj(i) & ~isnan(data(:,6)); % Get the data for that subject excluding non-responses
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
    
    params.nSims = 1000; % Set # of simulations per value difference 
    if type==0 % attribute-wise
        params.possDiffs = unique([ddmData.DA'; ddmData.DT']','rows'); % Set value differences experienced by participant
        params.type=0;
    else % option-wise
        params.possDiffs = unique([ddmData.AL'; ddmData.AR'; ddmData.TL'; ddmData.TR']','rows'); %Set value combos experienced
        params.type=1;
    end

    params.bias = 0; % starting bias. milosavljevic 2010 was -0.0002.
    params.rvsNoise = .1; 
    RTBins=quants;
    respOpts = [-1, 1];

    params.d1 = ddm(i,1); % Amt
    params.d2 = ddm(i,2); % time
    params.ndt1 = ddm(i,3); % Amt
    params.ndt2 = ddm(i,4); % Time
    params.bound = ddm(i,5);
     
    out = mtDDM_TestParams(ddmData,params,RTBins,respOpts);
        
    %Make RT distributions
    origQuants(i,:)=quantile(data(ind,7),linspace(.025,.975,8));

    sizes=size(out.rt);
    sampleRts=reshape(out.rt,[sizes(1)*sizes(2),1])/100;
    sampleRts(sampleRts==0)=[];%get rid of simulated non-responses
    sampleQuants(i,:)=quantile(sampleRts,linspace(.025,.975,8));
    
    if type==0 %attribute-wise
        imps=find(params.possDiffs(:,1)<0,1,'last'); %point below which right>left amt
        pats=find(params.possDiffs(:,1)>0,1,'first'); %point above which left>right amt
    else %option-wise
        imps=find(params.possDiffs(:,1)<params.possDiffs(:,2),1,'last'); %point below which right>left amt
        pats=find(params.possDiffs(:,1)>params.possDiffs(:,2),1,'first'); %point above which left>right amt
    end
    
    imp=length(find(out.choice(:,1:imps)==1))+length(find(out.choice(:,pats:end)==-1));
    pat=length(find(out.choice(:,1:imps)==-1))+length(find(out.choice(:,pats:end)==1));
    estPat(i,1)=pat/(imp+pat);
    actPat(i,1)=mean(data(ind,6));   
end
simRTs=[origQuants sampleQuants];
simAcc=[actPat estPat];
if sample==1;
    if type==0; %attribute wise
        csvwrite('attSimRT.csv',simRTs)
        csvwrite('attSimAcc.csv',simAcc)
    else %option-wise
        csvwrite('optSimRT.csv',simRTs)
        csvwrite('optSimAcc.csv',simAcc)
    end
else
    if type==0; %attribute wise
        csvwrite('attSimRT_rep.csv',simRTs)
        csvwrite('attSimAcc_rep.csv',simAcc)
    else %option-wise
        csvwrite('optSimRT_rep.csv',simRTs)
        csvwrite('optSimAcc_rep.csv',simAcc)
    end
end