% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% loads AnimalTA tracking output and computes basic parameters
% produces various plots

%% load data
clear;


assayPath="/Path/To/AnimalTA/Project/Results/Detailed_data";

%load chamber parameters
tmp=readtable("/Path/To/AnimalTA/Project/Results/Spatial/Arenas_coordinates1.csv");

generalCSV=readtable("/Path/To/AnimalTA/Project/Results/Spatial/General.csv");


loadIndVelThr="/Path/To/AnimalTA/Project/Results/Results_by_ind.csv";

decFflip=0;


%% process data

assayLength=30*60*20; %number of frames

baseParams=table2struct(tmp);


assayData=struct();
c=1;
%load and save parameters
for k=1:numel({baseParams.Video})
    folderPath=strcat(assayPath,'/',baseParams(k).Video);
    filePattern = sprintf('%s/*.csv', folderPath);
    baseFileNames = dir(filePattern);
    baseFileNames=baseFileNames(~startsWith({baseFileNames.name},'.'));
    numFiles=size({baseFileNames.name},2);

    for i=1:numFiles
        animal=['Ind',num2str(i-1)];
        fullFilePath = fullfile(folderPath, baseFileNames(i).name);
        tmpdata = readtable(fullFilePath);
        
        assayData(c).Video=baseParams(k).Video;
        assayData(c).Individual=animal;
        assayData(c).Center_X=baseParams(k).Center_X;
        assayData(c).Center_Y=baseParams(k).Center_Y;
        assayData(c).Xmin=baseParams(k).Xmin;
        assayData(c).Xmax=baseParams(k).Xmax;
        assayData(c).Ymin=baseParams(k).Ymin;
        assayData(c).Ymax=baseParams(k).Ymax;
        
        assayData(c).Surface=baseParams(k).Surface;
        assayData(c).Perimeter=baseParams(k).Perimeter;

        
        

        assayData(c).time=tmpdata.Time;
        assayData(c).x=tmpdata.X;
        assayData(c).y=tmpdata.Y;
        assayData(c).vel=tmpdata.Speed;
    
        %check what type of assay - acv=1, agar=0
        if contains(assayData(c).Video,'acv','IgnoreCase',true)
            assayData(c).acv=1;
            assayData(c).agar=0;
        elseif contains(assayData(c).Video,'agar','IgnoreCase',true)
            assayData(c).acv=0;
            assayData(c).agar=1;
        end
    
       
        %species are coded by number: ere=1, mel=2
        if contains(assayData(c).Video,'ere')
            assayData(c).spec=1;
        elseif contains(assayData(c).Video,'mel')
            assayData(c).spec=2;
        elseif contains(assayData(c).Video,'csMai')
            assayData(c).spec=2;
        end
        %sex is also coded by number: male=1, female=0
        if contains(assayData(c).Video,'WTm', 'IgnoreCase', true)
            assayData(c).sex=1;
        elseif contains(assayData(c).Video,'WTf', 'IgnoreCase', true)
            assayData(c).sex=0;
        end
        c=c+1;
    end
end

% add other parameters such as moving threshold

tmp=readtable(loadIndVelThr);
 
for i=1:size(tmp,1)

    target=find(matches(cellstr(str2mat(assayData.Video)),tmp.Video(i)) & matches(cellstr(str2mat(assayData.Individual)),tmp.Individual(i)));
    assayData(target).Moving_threshold=tmp.Moving_threshold(i);
    %assayData(target).recording_rig=tmp.recording_rig(i);
end

%correct velocity for moving threshold
for i=1:size(assayData,2)
    assayData(i).vel=assayData(i).vel-assayData(i).Moving_threshold;
    assayData(i).vel(assayData(i).vel<=0)=0;
    assayData(i).meanVel=mean(assayData(i).vel, 'omitnan');
end

%normalize coordinate system
for i=1:size(assayData,2)

    tmpx=assayData(i).x(1:assayLength)-assayData(i).Center_X;
    tmpy=assayData(i).y(1:assayLength)-assayData(i).Center_Y;
    assayData(i).radius=assayData(i).Perimeter/(2*pi);
    assayData(i).xnorm=tmpx./assayData(i).radius;
    assayData(i).ynorm=tmpy./assayData(i).radius;
    
end


gen=table2struct(generalCSV);

current={};
c=1;
% add bead locations and correct
for i=1:numel(gen)

    if i==1
        current=gen(i).Video;
        beadPos(c).Video=current;
        beadPos(c).beadx=[];
        beadPos(c).beady=[];
        beadPos(c).rawx=[];
        beadPos(c).rawy=[];
    else
        if ~isequal(gen(i).Video, current)
            current=gen(i).Video;
            c=c+1;
            beadPos(c).Video=current;
            beadPos(c).beadx=[];
            beadPos(c).beady=[];
            beadPos(c).rawx=[];
            beadPos(c).rawy=[];
        end
    end

    beadPos(c).beadx=[beadPos(c).beadx (gen(i).Element_Center_X-assayData(c).Center_X)/assayData(c).radius];
    beadPos(c).beady=[beadPos(c).beady (gen(i).Element_Center_Y-assayData(c).Center_Y)/assayData(c).radius];

    beadPos(c).rawx=[beadPos(c).rawx gen(i).Element_Center_X];
    beadPos(c).rawy=[beadPos(c).rawy gen(i).Element_Center_Y];
    
end



%% rough x,y plotting of all tracks in all assays



titles={'agar (n=5x5)','acv (n=5x5)'};
figure
for i=1:numel(groups)


    subplot(1,2,i)
    for j=groups(i).groups
        hold on
        scatter(assayData(j).xnorm,assayData(j).ynorm,6,'k','filled','MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05)

    end
    title(titles{i})
    xlim([-1 1])
    ylim([-1 1])
    axis off

    if i==2
        beadFind=find(contains({beadPos.Video},'acv','IgnoreCase',true));

    else
        beadFind=find(contains({beadPos.Video},'agar','IgnoreCase',true));
    end

    for j=1:length(beadFind)
        hold on
       % scatter(beadPos(j).beadx,beadPos(j).beady,'filled') %add beads to
       % plot if needed


    end
   
    
end



%% plot individual assay
%select assay
tmp=find(contains({assayData.Individual},'Ind0')); 
toplot=tmp([assayData(tmp).acv]==1);

%alternative
%toplot=find(contains({assayData.Video}, '012725_ere_1WTm_1_8do_Acv1DeadFemaleDotGridPos7_concat.avi'));

for i=toplot
    figure

    hold on

    hits=find(contains({assayData.Video},assayData(i).Video));
    
    fid=hits(hits~=i);
    

    scatter(assayData(i).xnorm,assayData(i).ynorm,6,'k','filled','MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05)
    scatter(assayData(fid).xnorm,assayData(fid).ynorm,6,'r','filled','MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05)
    
    
end
xlim([-1 1])
ylim([-1 1])
axis off

%% plot individual male trajectories and female position of 1m 1decapitated female assays

% Define the number of colors for each transition segment
nColors = 256*2;
nSegment = round(nColors / 2);  % Half for white-to-green, half for green-to-purple

% Transition from white to dark green
whiteToGreen = [linspace(1, 0, nSegment)', linspace(1, 0.5, nSegment)', linspace(1, 0, nSegment)'];

% Transition from dark green to purple
greenToPurple = [linspace(0, 0.5, nSegment)', linspace(0.5, 0, nSegment)', linspace(0, 0.5, nSegment)'];


% Combine the segments to form the full colormap
cmap = [whiteToGreen; greenToPurple];

%extract assays
tmp=find(contains({assayData.Individual},'Ind0'));
toplot=tmp([assayData(tmp).acv]==1);

 i=toplot(2); %define assay to plot (currently #2)
figure

hold on

hits=find(contains({assayData.Video},assayData(i).Video));

fid=hits(hits~=i);


 
% colored line

t = linspace(0,20,36000)'; % time in minutes

x=assayData(i).x(1:36000);
y=assayData(i).y(1:36000);

%add bead position
scatter(beadPos(2).rawx,beadPos(2).rawy,[],'r','filled')

plot(x, y,'k')


axis equal off

%% plot distance to closest bead and female over time
%1 - to closest bead
%2 - to female

%extract individual assays from stored list
tmp=find(contains({assayData.Individual},'Ind0'));

%define position of female in grid
fpos=[];

c=1;
for h=0
    figure()
    all_valuesF=[];
    all_valuesB=[];
    if h==1
        toplot=tmp([assayData(tmp).acv]==0);
    elseif h==0
        toplot=tmp([assayData(tmp).acv]==1);
    end
    for g=1:numel(toplot)
        i=toplot(g);

        focV=find(contains({assayData.Video},assayData(i).Video,'IgnoreCase',true));
        f_is=find(contains({assayData(focV).Individual},'Ind1'));
        fem=focV(f_is);

        
        %plot exactly 20 minutes from the beginning of the assay (30 fps
        %video recording)
        x=assayData(i).x(1:36000);
        y=assayData(i).y(1:36000);
        xf=assayData(fem).x(1:36000);
        yf=assayData(fem).y(1:36000);

        vid=find(contains({beadPos.Video},assayData(g).Video,'IgnoreCase',true));

        bx=beadPos(vid).rawx;
        by=beadPos(vid).rawy;
        
        P=[bx; by]'; %bead position
        
        

        assayData(i).pdistBeadsF=nan(numel(beadPos(vid).beadx),j);
        for j=1:numel(x)
        
            p=[x(j), y(j)];
            f=[xf(j),yf(j)];

            % Compute distances
            distances = sqrt(sum((P - p).^2, 2));
            mfdist=sqrt(sum((f - p).^2));

            

            assayData(i).pdistBeadsF(:,j)=distances;
            assayData(i).mfdist(:,j)=mfdist;
        
        end
        fwhere=fpos(c);
        c=c+1;
        female=assayData(i).mfdist;
        
        beads=assayData(i).pdistBeadsF;
        beads(fwhere,:)=[];
        
        
        
        %translate units to cm
        
        unit_diameter=assayData(i).Perimeter/pi;
        
        mm_diameter=190;
        ratio=unit_diameter/mm_diameter;
        
        xadj=0:35999;
        if ~isempty(assayData(i).firstApproachF)
            xadj=xadj-assayData(i).firstApproachF;
        
            from=assayData(i).firstApproachF-((30*120)+1);
            to=assayData(i).firstApproachF+((30*120)-1);
            
            lensT=30*120*2+1;
            sdistB=min(beads);
            if from<1
                from=1;
                tmp1=female(from:to);
                
                tmp2=sdistB(from:to);
                lens1=numel(tmp1);
                lensm=lensT-lens1;
                nn=nan(1,lensm);
                insertF=[nn tmp1];
                insertB=[nn tmp2];
            else
                insertF=female(from:to);
                insertB=sdistB(from:to);

            end

            all_valuesF=[all_valuesF; insertF];
            all_valuesB=[all_valuesB; insertB];
        subplot(2,1,1)
        hold on
        plot((-30*120):(30*120),insertB*ratio,'k')
        saver(g).bead=sdistB*ratio;
        title('distance to closest bead')
        ylabel('distance [mm]')
        xlim([-30*120,30*120])
        ylim([0 60])
        yline(3,'--r')
        
        subplot(2,1,2)
        hold on
        plot((-30*120):(30*120),insertF*ratio,'k')
        title('distance to female')
        
        xticks(0:30*60*2:36000)
        xticklabels({'0','2','4','6','8','10','12','14','16','18','20'})
        xlim([-30*120,30*120])
        ylabel('distance [mm]')
        yline(3,'--r')
        ylim([0 120])
        xlabel('time [min]')
        
        end
        

    end

    subplot(2,1,2)
    hold on
    plot((-30*120):(30*120),mean(all_valuesF,'omitnan')*ratio,'k','linewidth',2)
    subplot(2,1,1)
    hold on
    plot((-30*120):(30*120),mean(all_valuesB,'omitnan')*ratio,'k','linewidth',2)

end

%% plot heatmap of positions in assay chamber

% Define the number of colors for each transition segment
nColors = 256*2;
nSegment = round(nColors / 2);  % Half for white-to-green, half for green-to-purple

% Transition from white to dark green
whiteToGreen = [linspace(1, 0, nSegment)', linspace(1, 0.5, nSegment)', linspace(1, 0, nSegment)'];

% Transition from dark green to purple
greenToPurple = [linspace(0, 0.5, nSegment)', linspace(0.5, 0, nSegment)', linspace(0, 0.5, nSegment)'];


% Combine the segments to form the full colormap
cmap = [whiteToGreen; greenToPurple];

% iterate over all assays in predefined group and calculate relative
% position

for i=1:numel(groups)
    figure()
    xall=vertcat(assayData(groups(i).groups).xnorm);
    yall=vertcat(assayData(groups(i).groups).ynorm);
    
    % Define bin edges
    xEdges = linspace(min(xall), max(xall), 25); % 50 bins along x-axis
    yEdges = linspace(min(yall), max(yall), 25); % 50 bins along y-axis
    
    % Bin the data using histcounts2
    [counts, xEdges, yEdges] = histcounts2(xall, yall, xEdges, yEdges);
    
    % Total number of points (for relative occupancy)
    totalPoints = sum(counts(:));
    
    % Normalize to relative occupancy
    relativeOccupancy = counts / totalPoints;
    max(max(relativeOccupancy))
    % Visualize normalized heatmap
    imagesc(xEdges, yEdges, relativeOccupancy');
    saver(i).relOcc=relativeOccupancy';
    axis xy; % Correct axis orientation
    colormap(cmap); % Choose a colormap
    colorbar; % Add a color bar for reference
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Heatmap of XY Coordinates');
    xlim([-1 1])
    ylim([-1 1])
    %caxis([0 0.02])
end
%% infer heading, appraoch, visits

% calculates heading angle for each individual
% calculates the pdistance between the animal at each
% point to all beads in the assay chamber
timeThreshold=30; %in frames
distThreshold=0.05;
for i=1:numel(assayData)
    
    %only for 1m1decapf assays
    if decFflip==1
        if contains({assayData(i).Individual},'Ind1'), continue;end %skips decap F for stats
    end

    x=assayData(i).xnorm;
    y=assayData(i).ynorm;

    whichBead=find(matches({beadPos.Video},assayData(i).Video));

    if numel(beadPos(whichBead).beadx)>13
        max=numel(beadPos(whichBead).beadx);
        while max>13
            beadPos(whichBead).beadx(max)=[];
            beadPos(whichBead).beady(max)=[];
            beadPos(whichBead).rawx(max)=[];
            beadPos(whichBead).rawy(max)=[];


            max=max-1;
        end
    end

    xBeads=beadPos(whichBead).beadx;
    yBeads=beadPos(whichBead).beady;

    distances=nan(length(x),numel(xBeads));
    numBeads=0; %close to how many beads during assay?
    numCloseFrames=0;



    for j=1:numel(xBeads)
        
        x0=xBeads(j);
        y0=yBeads(j);

        %euclidean distances
        dist=sqrt((x - x0).^2 + (y - y0).^2);
        distances(:,j) = dist;

        closeFrames=sum(dist<distThreshold);
        if closeFrames>timeThreshold
            numBeads=numBeads+1;
        end
        numCloseFrames=numCloseFrames+closeFrames;
        
    end
    assayData(i).distToBeads=distances;
    assayData(i).beadsVisited=numBeads;
    assayData(i).framesAtBead=numCloseFrames;


    %identify female position (only use for decapitated female assays -> use decFflip switch)
    if decFflip==1
        hits=find(contains({assayData.Video},assayData(i).Video));
        
        fid=hits(hits~=i);
    
        xf=assayData(fid).xnorm;
        yf=assayData(fid).ynorm;
    
        
        fdist=sqrt((x - xf).^2 + (y - yf).^2);
        fclose=sum(fdist<distThreshold);
    
        fframes=find(fdist<distThreshold);
    
        assayData(i).distToF=fdist;
        assayData(i).numCloseFramesF=fclose;
        if fframes
            assayData(i).firstApproachF=fframes(1);
        else
            assayData(i).firstApproachF=[];
        end

    end


    heading = zeros(length(x)-1);

    % Compute differences between consecutive coordinates
    dx=diff(x);
    dy=diff(y);

    % Compute heading in radians using atan2
    theta_rad = atan2(dy, dx);

    assayData(i).heading=theta_rad;


end

% extracts number of visits for each individual (including repeat visits) and the respective
% lengths of said visits
for i=1:numel(assayData)

    threshold=distThreshold; % threshold distance to call individual 'at bead'
    tolerance=15; %allow for half a second tolerance in consecutive time at bead

    A=assayData(i).distToBeads;

    %loop bead by bead
    visitsLens=[];
    visitsBead=[];
    approach=[];
    for j=1:size(assayData(i).distToBeads,2)
        
        a0=A(:,j); %focal bead

        %create logical mask where true means a0 is below threshold
        mask = a0<threshold;

        %tolerance vector
        se = ones(tolerance,1);

        mask_filled = imclose(mask, se);

        %identify contiguous regions in the filled mask
        [labels, numSegments] = bwlabel(mask_filled);


        segments = cell(numSegments, 1);
        for k = 1:numSegments
            indices = find(labels == k);
            segments{k} = indices;  % or a0(indices) to extract the actual data points
            %fprintf('Segment %d: Start index = %d, End index = %d\n', j, indices(1), indices(end));
            if length(indices)>=timeThreshold
                visitLens=length(indices);
                visitsLens=[visitsLens visitLens];
                visitsBead=[visitsBead j];
                approach=[approach indices(1)];
            end

        end
    end

    assayData(i).visitLens=visitsLens;
    assayData(i).visitBeadID=visitsBead;
    assayData(i).numVisits=length(visitsBead);


    [sorted_approach, sort_idx]=sort(approach);
    assayData(i).BeadApproaches=sorted_approach;
    assayData(i).BeadIdsSorted=visitsBead(sort_idx);


end



%% plot x,y and heading (to check performance)

for i=1:numel(assayData)
    
    %only for 1m1decapf assays
    if contains({assayData(i).Individual},'Ind1'), continue;end %skips decap F for stats




    if assayData(i).acv ==1
        col='r';
    else
        col='b';
    end
    heading=assayData(i).heading;
    x=assayData(i).xnorm;
    y=assayData(i).ynorm;

    arrowLength = 0.05;

    sampleRate = 10;
    figure;
    quiver(x(1:sampleRate:end-1), y(1:sampleRate:end-1), arrowLength*cos(heading(1:sampleRate:end)), ...
    arrowLength*sin(heading(1:sampleRate:end)), 0, 'MaxHeadSize', 0.5, 'Color', col);
    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Coordinates with Heading Arrows (sampled every 10 frames)');
    hold on
    scatter(beadPos(i).beadx,beadPos(i).beady,50,'k','filled')
    xlim([-1 1])
    ylim([-1 1])
end

%% plot bead visitation stats

figure
subplot(2,2,1)
for i=1:numel(groups)
    
    
    y=[assayData(groups(i).groups).beadsVisited];
    x=i*ones(size(y));
    saver(i).beadsVisited=y';

    hold on
    scatter(x,y,'filled','jitter','on','jitterAmount',0.1)
    plot([i-.15 i+.15],[mean(y) mean(y)],'k','LineWidth',2)
    
    xlim([0.75 2.25])
    ylim([0 13])
    xticks([1 2])
    xticklabels({'agar','acv'})
    ylabel('Beads visited')

end

subplot(2,2,2)
for i=1:numel(groups)

    y=[assayData(groups(i).groups).numVisits];
    
    x=i*ones(size(y));
    saver(i).numVisits=y'/20;

    hold on
    scatter(x,y,'filled','jitter','on','jitterAmount',0.1)
    plot([i-.15 i+.15],[mean(y) mean(y)],'k','LineWidth',2)
    
    xlim([0.75 2.25])
    %ylim([0 13])
    xticks([1 2])
    xticklabels({'agar','acv'})
    ylabel('Total visits')


end

subplot(2,2,3)
for i=1:numel(groups)

    y=[assayData(groups(i).groups).framesAtBead];
    y=y/30;

    x=i*ones(size(y));
    hold on
    scatter(x,y,'filled','jitter','on','jitterAmount',0.1)
    plot([i-.15 i+.15],[mean(y) mean(y)],'k','LineWidth',2)
    
    xlim([0.75 2.25])
    %ylim([0 13])
    xticks([1 2])
    xticklabels({'agar','acv'})
    ylabel('Total time at beads [sec]')



end

%visit length
subplot(2,2,4)
for i=1:numel(groups)

    y=[assayData(groups(i).groups).visitLens];
    y(y<30)=[]; %remove visits below a second
    y=y/30;

    x=i*ones(size(y));
    hold on
    scatter(x,y,'filled','jitter','on','jitterAmount',0.1)
    plot([i-.15 i+.15],[mean(y) mean(y)],'k','LineWidth',2)
    
    xlim([0.75 2.25])
    xticks([1 2])
    xticklabels({'agar','acv'})
    ylabel('Visit length [sec]')
end

%% plot approaches per bead
figure
titles={'mean visits per bead agar','mean visits per bead acv','mean total time at bead agar','mean total time at bead acv'};


% Define the number of colors for each transition segment
nColors = 256*2;
nSegment = round(nColors / 2);  % Half for white-to-green, half for green-to-purple

% Transition from white to dark green
whiteToGreen = [linspace(0.9, 0, nSegment)', linspace(0.9, 0.5, nSegment)', linspace(0.9, 0, nSegment)'];

% Transition from dark green to purple
greenToPurple = [linspace(0, 0.5, nSegment)', linspace(0.5, 0, nSegment)', linspace(0, 0.5, nSegment)'];


% Combine the segments to form the full colormap
cmap = [whiteToGreen; greenToPurple];



for i=1:numel(groups)

    
    ind=groups(i).groups;
    counts=zeros(13,length(ind));
    time=counts;
    switches=zeros(size(ind));
    for j=1:length(ind)

        counts_loop=zeros(13,1);
        time_loop=counts_loop;

        for k=1:13
            counts_loop(k)=sum(assayData(ind(j)).visitBeadID==k);
            time_loop(k)=sum(assayData(ind(j)).visitLens(assayData(ind(j)).visitBeadID==k));
        end

        counts(:,j)=counts_loop;
        time(:,j)=time_loop;
        switches(j)=sum(diff(assayData(ind(j)).BeadIdsSorted)~=0);
    end
    
    meanCounts=mean(counts,2);
    meanTime=mean(time,2);

    groups(i).meanBeadVisits=meanCounts;
    groups(i).totalBeadVisits=counts;
    groups(i).meanTimeAtBead=meanTime;
    groups(i).timeByInd=time;
    groups(i).switches=switches;

    

    subplot(2,2,i)

    x=beadPos(1).beadx;
    y=beadPos(1).beady;
    
    scatter(x,y,6,'k')
    hold on
    sizeP=50*meanCounts;
    sizeP(sizeP==0)=nan();

    scatter(x,y,sizeP,round(meanCounts),'filled')
    colormap(cmap)
    colorbar
    caxis([0 8])
    xlim([-1 1])
    ylim([-1 1])
    axis off

    title(titles{i})
    subplot(2,2,i+2)

    x=beadPos(1).beadx;
    y=beadPos(1).beady;
    
    scatter(x,y,6,'k')
    hold on
    sizeP=meanTime/20;
    sizeP(sizeP==0)=nan();

    scatter(x,y,sizeP,round(meanTime/30),'filled')
    colormap(cmap)
    colorbar
    caxis([0 60])
    xlim([-1 1])
    ylim([-1 1])
    axis off

    title(titles{i+2})
end

figure
for i=1:numel(groups)

    y=groups(i).switches;
    x=i*ones(size(y));

    swm=mean(y);
    hold on
    scatter(x,y,[],"filled",'jitter','on')

    plot([i-0.25 i+0.25],[swm swm],'k','LineWidth',2)
end

xlim([0.5 2.5])
ylim([0 4*20])


%% plot approaches per individual
figure
titles0={'bead approaches agar', 'bead approaches acv'};
for i=1:numel(groups)
    subplot(1,2,i)
    hold on
    c=1;
    c2=c;
    ind=groups(i).groups;

    for j=1:length(ind)

        x=assayData(ind(j)).BeadApproaches;
        saver(ind(j)).beadApproaches=strjoin(string(x/30), ',');
        y=c2*ones(size(x));
        scatter(x,y,[],'k','filled', 'MarkerFaceAlpha', 0.45)

        if mod(c,5)==0
            yline(c2+1,'k')
            c2=c2+2;
            
        else
            c2=c2+1;
        end
        c=c+1;
    end


    ylim([0 30])
    yticks(1:6:c2-1)
    yticklabels(repmat({'male 1'}, 1, round(length(ind)/5)))
    xlim([0 assayLength])
    xticks(0:30*120:assayLength)
    xticklabels(num2cell(0:2:20))
    xlabel('time (min)')
    title(titles0{i})
end



%%

% extract bead approaches pre vs post female encounter
% extract male time spent at each bead vs female

for i=1:numel(assayData)


    %get female id
    str=assayData(i).Video;
    tokens = regexp(str, 'GridPos(\d+)', 'tokens');

    if ~isempty(tokens)
        numberAfterGridPos = str2double(tokens{1}{1});
    else
        numberAfterGridPos = NaN;  % or handle the case where 'GridPos' is not found
    end

    assayData(i).femalePos=numberAfterGridPos;


    timePerBead=zeros(13,1);
    timePerBeadpre=timePerBead;
    timePerBeadpost=timePerBead;
    for j=1:13 %number beads

        if j==assayData(i).femalePos
            
            timePerBead(j)=nan();
            timePerBeadpre(j)=nan();
            timePerBeadpost(j)=nan();
            continue;
        end

        bint=find(assayData(i).visitBeadID==j);
        ttime=sum(assayData(i).visitLens(bint));
        timePerBead(j)=ttime;

        if assayData(i).firstApproachF

            pre=find(assayData(i).BeadApproaches<assayData(i).firstApproachF);
            post=find(assayData(i).BeadApproaches>assayData(i).firstApproachF);

            bintpre=find(assayData(i).visitBeadID(pre)==j);
            bintpost=find(assayData(i).visitBeadID(post)==j);

            ttimepre=sum(assayData(i).visitLens(bintpre));
            ttimepost=sum(assayData(i).visitLens(bintpost));
            timePerBeadpre(j)=ttimepre;
            timePerBeadpost(j)=ttimepost;

        end


    end
   

    assayData(i).bead_int_totalTime=timePerBead;
    assayData(i).bead_int_totalTimepre=timePerBeadpre;
    assayData(i).bead_int_totalTimepost=timePerBeadpost;
end

%% plot bead approaches pre vs post female encounter

tmp=find(contains({assayData.Individual},'Ind0'));
toplot=tmp([assayData(tmp).acv]==1);


a=[assayData(toplot).bead_int_totalTime]/assayLength;



div=zeros(size([assayData(toplot).numCloseFramesF]));
div2=repmat(assayLength,size([assayData(toplot).numCloseFramesF]));

c=1;
for i=toplot
    if assayData(i).firstApproachF
        div(c)=(assayLength-[assayData(i).firstApproachF])+1;
    end
    c=c+1;
end
div2(1:length([assayData(toplot).firstApproachF]))=[assayData(toplot).firstApproachF]-1;
pre=[assayData(toplot).bead_int_totalTimepre]./div2;

post=[assayData(toplot).bead_int_totalTimepost]./div;

figure
hold on

scatter(ones(size(mean(a,'omitnan'))),mean(a,'omitnan'),'filled','jitter','on','jitteramount',0.1)
plot([0.75 1.25], [mean(a(:),'omitnan') mean(a(:),'omitnan')],'k','linewidth',2)
scatter(4*ones(size([assayData(toplot).numCloseFramesF])),[assayData(toplot).numCloseFramesF]./div,'filled','jitter','on','jitteramount',0.1)
plot([3.75 4.25], [mean([assayData(toplot).numCloseFramesF]./div,'omitnan') mean([assayData(toplot).numCloseFramesF]./div,'omitnan')],'k','linewidth',2)

scatter(2*ones(size(mean(pre,'omitnan'))),mean(pre,'omitnan'),'filled','jitter','on','jitteramount',0.1)
plot([1.75 2.25], [mean(pre(:),'omitnan') mean(pre(:),'omitnan')],'k','linewidth',2)

scatter(3*ones(size(mean(post,'omitnan'))),mean(post,'omitnan'),'filled','jitter','on','jitteramount',0.1)
plot([2.75 3.25], [mean(post(:),'omitnan') mean(post(:),'omitnan')],'k','linewidth',2)



xlim([0.5 4.5])
xticks(1:4)
xticklabels({'per bead','pre f encounter','post f encounter','f post encounter'})
ylabel('average fraction of time')
