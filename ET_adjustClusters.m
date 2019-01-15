function ET_adjustClusters(sample)
% Original version written by Dianna Amasino 2016
% Updated by Dianna Amasino 2019 based on edits made by Khoi Vo, 2016
% Huettel lab, Duke University

% This function takes in eye tracking data classified by Tobii as fixations 
% and extracted from the choice period of trials and adjusts it so that the 
% clusters relating to each AOI are properly centered (accounting for bad 
% calibration or calibration drift over the course of the experiment. 

% It does this by:
% 1) Fitting clusters of the data.
% 2) Identifying clusters at the 4 corners (which excludes clusters left 
% over from fixation).
% 3) Moving the clusters such that the center of mass of the cluster is 
% centered in the AOI.

% INPUT
% sample - 1 (primary) or 2 (replication) data
% col 1: subject number (matching behavioral data)
% col 2: Trial number
% col 3: x-position (pixels on screen)
% col 4: y-position (pixels on screen)
% col 5: fixation counter

% Screen dimensions
% width_x = 1280;
% height_y = 1024; 

ET_adj=[];

%Note that in psychtoolbox, (0,0) is the TOP left of the screen
aoi.TL = [65 416 350 350]; % Top left AOI
aoi.BL = [65 1016 350 350]; % Bottom left AOI
aoi.TR = [865 416 350 350]; % Top right AOI
aoi.BR = [865 1016 350 350]; % Bottom right AOI
aois=fieldnames(aoi);

xpos.TL= 240; % Middle of left aoi, x pos
xpos.TR= 1040; % Middle of right aoi, x pos
xpos.BL= 240; % Middle of left aoi, x pos
xpos.BR= 1040; % Middle of right aoi, x pos
ypos.TL= 241; % Middle of top aoi, y pos
ypos.BL= 841; % Middle of bottom aoi, y pos
ypos.TR= 241; % Middle of top aoi, y pos
ypos.BR= 841; % Middle of bottom aoi, y pos
 
% maximum cluster to extract, can edit
maxcluster = 14;

%Load data
dataPath=pwd; %adapt to your location
cd(dataPath)
if sample ==1 %Primary sample
    load('amasinoETAl_ET.csv') 
    data=amasinoETAl_ET;
    subj=[1:5,7:17,19,21:56,58,60:63,65,68,70:73,75:91,93:101,103:117]; %accounting for ET exclusions
else % replication sample
    load('amasinoETAl_ET_rep.csv')
    data=amasinoETAl_ET_rep;
    subj=[1:3,5:9,11:17,19,21:24,26:38,40,41,43:49,51:55,58,60,62:71,73:81,83:94,96:100];
end

%Loop over subjects
for i = 1:length(subj)
    adjET=[];
    adjET1=[];
    adjET2=[];
    sub = data(:,1)==subj(i); %Find subject-specific data points
    grandmean=mean(data(sub,3:4)); %Find center of all eye data points
    
    % Cluster data according to all data (in case of slighty off calibration)
    clusters=clusterdata(data(sub,3:4),'linkage','centroid','savememory','on','maxclust',maxcluster);
    
    clusterMat=[clusters data(sub,:)]; %Concatenate clusters and data
    
    %Summarize info about clusters including cluster number, number of data
    %points in each cluster, and mean cluster x and y position
    means=grpstats(clusterMat,clusters);
    clusterInfo=[(1:14)',accumarray(clusters,1),means(:,4:5)];
        
    %Assign clusters to each quadrant of the screen relative to the grandmean as center
    %Clusters should be at least 100 pixels horizontally away from the center to avoid the fixation cluster
    ind.TL = clusterInfo(clusterInfo(:,3)<grandmean(1)-100 & clusterInfo(:,4)<grandmean(2),1);
    ind.BL = clusterInfo(clusterInfo(:,3)<grandmean(1)-100 & clusterInfo(:,4)>grandmean(2),1);
    ind.TR = clusterInfo(clusterInfo(:,3)>grandmean(1)+100 & clusterInfo(:,4)<grandmean(2),1);
    ind.BR = clusterInfo(clusterInfo(:,3)>grandmean(1)+100 & clusterInfo(:,4)>grandmean(2),1);
    
    %Loop until there is 1 and only 1 cluster (largest size) in quadrant for the 4 AOIs
    while length(ind.TL)~=1 | length(ind.BL)~=1 | length(ind.TR)~=1 | length(ind.BR)~=1
        for j=1:length(aois)
            if ~isempty(ind.(aois{j})) %if at least 1 cluster falls within quadrant
                [M,I]=max(clusterInfo(ind.(aois{j}),2)); %Find largest sized cluster
                % Find left-most cluster
            end
            ind.(aois{j})=ind.(aois{j})(I);
        end    
    end
    %Do the same as above, but separating the 1st (trials 1-70) and 2nd (trials 71-141) 
    % runs to separately correct for drift across runs if need be. 
    %Make clusters
    sub1 = data(:,1)==subj(i) & data(:,2)<71;
    sub2 = data(:,1)==subj(i) & data(:,2)>70;
    clusters1=clusterdata(data(sub1,3:4),'linkage','centroid','savememory','on','maxclust',maxcluster);
    clusters2=clusterdata(data(sub2,3:4),'linkage','centroid','savememory','on','maxclust',maxcluster);
    
    %concatenate clusters 
    clusterMat1=[clusters1 data(sub1,:)];
    clusterMat2=[clusters2 data(sub2,:)];
    
    %Summarize cluster info
    means=grpstats(clusterMat1,clusters1);
    clusterInfo1=[(1:14)',accumarray(clusters1,1),means(:,4:5)];
    means=grpstats(clusterMat2,clusters2);
    clusterInfo2=[(1:14)',accumarray(clusters2,1),means(:,4:5)];
    
    %Assign clusters to each quadrant for 1st run. Here, add extra distance
    %in vertical direction as well (more for bottom b/c fixation is closer to top)
    ind1.TL = clusterInfo1(clusterInfo1(:,3)<grandmean(1)-100 & clusterInfo1(:,4)<grandmean(2)-50,1);
    ind1.BL = clusterInfo1(clusterInfo1(:,3)<grandmean(1)-100 & clusterInfo1(:,4)>grandmean(2)+100,1);
    ind1.TR = clusterInfo1(clusterInfo1(:,3)>grandmean(1)+100 & clusterInfo1(:,4)<grandmean(2)-50,1);
    ind1.BR = clusterInfo1(clusterInfo1(:,3)>grandmean(1)+100 & clusterInfo1(:,4)>grandmean(2)+100,1);
    
    %Loop until there is 1 and only 1 cluster (largest size) in quadrant for the 4 AOIs
    while length(ind1.TL)~=1 | length(ind1.BL)~=1 | length(ind1.TR)~=1 | length(ind1.BR)~=1
        for j=1:length(aois)
            if ~isempty(ind1.(aois{j})) %if at least 1 cluster falls within quadrant
                [M,I]=max(clusterInfo1(ind1.(aois{j}),2)); %Find largest sized cluster
                % Find left-most cluster
            end
            ind1.(aois{j})=ind1.(aois{j})(I);
        end    
    end
    clustInd1=[ind1.TL ind1.BL ind1.TR ind1.BR]; % concatenate all clusters for use below
    
    %Assign clusters to each quadrant for 2nd run. Here, add extra distance
    %in vertical direction as well (more for bottom b/c fixation is closer to top) 
    ind2.TL = clusterInfo2(clusterInfo2(:,3)<grandmean(1)-100 & clusterInfo2(:,4)<grandmean(2)-50,1);
    ind2.BL = clusterInfo2(clusterInfo2(:,3)<grandmean(1)-100 & clusterInfo2(:,4)>grandmean(2)+100,1);
    ind2.TR = clusterInfo2(clusterInfo2(:,3)>grandmean(1)+100 & clusterInfo2(:,4)<grandmean(2)-50,1);
    ind2.BR = clusterInfo2(clusterInfo2(:,3)>grandmean(1)+100 & clusterInfo2(:,4)>grandmean(2)+100,1);
    
    %Loop until there is 1 and only 1 cluster (largest size) in quadrant for the 4 AOIs
    while length(ind2.TL)~=1 | length(ind2.BL)~=1 | length(ind2.TR)~=1 | length(ind2.BR)~=1
        for j=1:length(aois)
            if ~isempty(ind2.(aois{j})) %if at least 1 cluster falls within quadrant
                [M,I]=max(clusterInfo2(ind2.(aois{j}),2)); %Find largest sized cluster
                % Find left-most cluster
            end
            ind2.(aois{j})=ind2.(aois{j})(I);
        end    
    end
    clustInd2=[ind2.TL ind2.BL ind2.TR ind2.BR]; % concatenate all clusters for use below
    
    % Adjust all gaze points within the AOI (350x350 pixels) such that the 
    % appropriate cluster is centered in the AOI.
    % If from run 1 to run 2 clusters are >40 pixels apart and clusters are
    % more than 50 gaze points, adjust separately, otherwise all at once.
   if sum(sum(abs(clusterInfo1(clustInd1,3:4)-clusterInfo2(clustInd2,3:4))>40)) ...
      && sum(clusterInfo1(clustInd1,2)>50) && sum(clusterInfo2(clustInd2,2)>50)
       
       for k=1:length(aois) %for each AOI  
           % Run 1: Find the points within the AOI (+/-175 pixels from center of cluster)
           adj1.(aois{k}) = find(clusterMat1(:,4)>clusterInfo1(ind1.(aois{k}),3)-175 & clusterMat1(:,4)<clusterInfo1(ind1.(aois{k}),3)+ 175 ...
           & clusterMat1(:,5)>clusterInfo1(ind1.(aois{k}),4)-175 & clusterMat1(:,5)<clusterInfo1(ind1.(aois{k}),4)+ 175);   
           
           % Adjust the x and y position such that cluster center is in middle of AOI
           adjET1(adj1.(aois{k}),1:6) = [clusterMat1(adj1.(aois{k}),2:3) ...
           clusterMat1(adj1.(aois{k}),4)+(xpos.(aois{k})-clusterInfo1(ind1.(aois{k}),3)) ...
           clusterMat1(adj1.(aois{k}),5)+(ypos.(aois{k})-clusterInfo1(ind1.(aois{k}),4)) ...
           clusterMat1(adj1.(aois{k}),6) repmat(k,length(adj1.(aois{k})),1)];
       
           % Run 2: Find the points within the AOI (+/-175 pixels from center of cluster)
           adj2.(aois{k}) = find(clusterMat2(:,4)>clusterInfo2(ind2.(aois{k}),3)-175 & clusterMat2(:,4)<clusterInfo2(ind2.(aois{k}),3)+ 175 ...
           & clusterMat2(:,5)>clusterInfo2(ind2.(aois{k}),4)-175 & clusterMat2(:,5)<clusterInfo2(ind2.(aois{k}),4)+ 175);   
           
           % Adjust the x and y position such that cluster center is in middle of AOI
           adjET2(adj2.(aois{k}),1:6) = [clusterMat2(adj2.(aois{k}),2:3) ...
           clusterMat2(adj2.(aois{k}),4)+(xpos.(aois{k})-clusterInfo2(ind2.(aois{k}),3)) ...
           clusterMat2(adj2.(aois{k}),5)+(ypos.(aois{k})-clusterInfo2(ind2.(aois{k}),4)) ...
           clusterMat2(adj2.(aois{k}),6) repmat(k,length(adj2.(aois{k})),1)];
       end
       adjET=[adjET1;adjET2]; % combine first and second runs 
        
   else %adjust together   
       for k=1:length(aois) %for each AOI          
           % Find the points within the AOI (+/-175 pixels from center of cluster)
           adj.(aois{k}) = find(clusterMat(:,4)>clusterInfo(ind.(aois{k}),3)-175 & clusterMat(:,4)<clusterInfo(ind.(aois{k}),3)+ 175 ...
           & clusterMat(:,5)>clusterInfo(ind.(aois{k}),4)-175 & clusterMat(:,5)<clusterInfo(ind.(aois{k}),4)+ 175);
           
           % Adjust the x and y position such that cluster center is in middle of AOI
           adjET(adj.(aois{k}),1:6) = [clusterMat(adj.(aois{k}),2:3) ...
           clusterMat(adj.(aois{k}),4)+(xpos.(aois{k})-clusterInfo(ind.(aois{k}),3)) ...
           clusterMat(adj.(aois{k}),5)+(ypos.(aois{k})-clusterInfo(ind.(aois{k}),4)) ...
           clusterMat(adj.(aois{k}),6) repmat(k,length(adj.(aois{k})),1)];
       end  
   end
   adjET(adjET(:,1)==0,:)=[]; %Get rid of points not in AOIs
   ET_adj=[ET_adj; adjET]; %Add this subject to all subjects
end

%Save output
if sample==1
    csvwrite('ET_adj.csv',ET_adj)
else
    csvwrite('ET_adj_rep.csv',ET_adj)
end