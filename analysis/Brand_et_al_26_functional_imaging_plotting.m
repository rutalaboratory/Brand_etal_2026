% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% loads and plots 2-photon calcium imaging data and tracking data saved as csv
% tracking data produced with Brand_et_al_26_functional_imaging_postProcessing.m
%
% produces various plots
%
% correct data loading requires the following naming scheme of files
% Imaging data: species_housing_stimulusParadigm_neuropil.csv
% Tracking data: species_housing_stimulusParadigm_neuropil_tracking.csv
%
% script loads all replicate experiments of a given experimental condition (as indicated by the name) and groups them


%% load data
clear;
start_path = fullfile('/Path/To/Data');
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end


filePattern = sprintf('%s/*.csv', topLevelFolder);
baseFileNames = dir(filePattern);
numberOfMatFiles = length(baseFileNames);
    
spec = {'species'};
housing = {'housing'};
chem={'stimulusParadigm'};
neuropils = {'neuropil'};


groups = strcat({'group'},spec,{'_'},housing, {'_'},chem,{'_'}, neuropils);
numGroups=length(groups);

patterns = strcat(spec,{'.*'},housing, {'.*'},chem,{'.*'},neuropils);

data=struct();

for pat = 1:length(patterns)
    maxLens=0; %save size of least number of frames recorded for interpolation
    maxTrack=0;

    thisPattern=patterns{pat};
    %fprintf('checking for %s\n',thisPattern)
    
    hits={};
    %get matches to pattern - these are assays within the same group

    for f = 1 : numberOfMatFiles
        if regexp(baseFileNames(f).name, thisPattern)
            %fprintf('\tfound %s\n',baseFileNames(f).name)
            hits{length(hits)+1} = baseFileNames(f).name;
        end
    end
        
    numHits=length(hits);
        
    for hit = 1:numHits
        fullFile = fullfile(topLevelFolder, hits{hit});
        fprintf('loading file %s\n',hits{hit})
        
        in=readtable(fullFile);
        con=table2array(in);
        hitref=sprintf('map%d',hit);
        fileRep.(groups{pat}).(hitref)=hits{hit};
        data.(groups{pat})(hit).fileName2p=hits{hit};
        lens=size(con,1);
        if lens>maxLens
            maxLens=lens;
            
        end

        raw.(groups{pat}).(hitref)=con;
        data.(groups{pat})(hit).f_raw=con(:,2);
        data.(groups{pat})(hit).fsize_raw=numel(con(:,2));

        %load tracking file
        id=split(hits{hit},'_');
        patt="m" + (digitsPattern(1)|digitsPattern(2)) + "v" + (digitsPattern(1)|digitsPattern(2));
        individual=find(contains(id,patt));
        trackPattern = strcat(id{1},{'.*'},spec,{'.*'},housing, {'.*'},chem,{'.*'},id{individual}, {'.*'},'tracking');

        for f = 1 : numberOfMatFiles
            if regexp(baseFileNames(f).name, trackPattern{1})
                %fprintf('\tfound %s\n',baseFileNames(f).name)
                trackFile = baseFileNames(f).name;
        
            end
        end

        fullFile = fullfile(topLevelFolder, trackFile);
        fprintf('loading file %s\n',trackFile)
        
        inT=readtable(fullFile);
        track=table2array(inT);
        data.(groups{pat})(hit).track_size_raw=size(track,1);
        data.(groups{pat})(hit).track_bTime_raw=track(:,1);
        data.(groups{pat})(hit).track_femaleAngle_raw=track(:,2);
        data.(groups{pat})(hit).track_hDiff_raw=track(:,3);
        data.(groups{pat})(hit).track_TI_raw=track(:,4);
        data.(groups{pat})(hit).track_fidelity_raw=track(:,5);
        data.(groups{pat})(hit).track_vigor_raw=track(:,6);
        data.(groups{pat})(hit).track_maleVel_raw=track(:,7);
        data.(groups{pat})(hit).track_maleAngVel_raw=track(:,8);


        tracking.(groups{pat}).(hitref)=track;
        lensT=size(track,1);
        if lensT>maxTrack
            maxTrack=lensT;
        end

    end
end

%interpolate 2p data
for i = 1:numGroups
    
   %totals first 
   numAssays = numel(data.(groups{i}));
   assays = fieldnames(raw.(groups{i}));

   longest=max([data.(groups{i}).fsize_raw]);
   XA=(1:longest)/(longest/210);

   for j = 1:numAssays
       focSize=data.(groups{i})(j).fsize_raw;
       data.(groups{i})(j).f_time_interp=XA;

       if focSize<longest
             
            B=data.(groups{i})(j).f_raw;
            XB=(1:focSize)/(focSize/210);
            
            new_B=interp1(XB,B,XA);
            data.(groups{i})(j).f_interp=new_B;
            

            interp_raw.(groups{i}).(assays{j})=[XA;new_B]';
            


       else
            data.(groups{i})(j).f_interp=data.(groups{i})(j).f_raw';
            interp_raw.(groups{i}).(assays{j})=[XA', raw.(groups{i}).(assays{j})(:,2)];
            
       end

       %interpolate tracking index to match number of frames recorded in 2p
 

       B=data.(groups{i})(j).track_TI_raw;
       XB=data.(groups{i})(j).track_bTime_raw;
       new_B=interp1(XB,B,XA);
        data.(groups{i})(j).track_TI_interp=new_B;
       

       %interpolate vis stim
       B2=data.(groups{i})(j).track_femaleAngle_raw;
       new_B2=interp1(XB,B2,XA);
       
       data.(groups{i})(j).track_femaleAngle_interp=new_B2;

       trackInterp.(groups{i}).(assays{j})=[XA; new_B; new_B2]';

       %calculate df/f

       fbase=mean(interp_raw.(groups{i}).(assays{j})(30:180,2)); %baseline 30:180
       fdelta=interp_raw.(groups{i}).(assays{j})(:,2)-fbase;
       fof=fdelta/fbase;
       dff.(groups{i}).(assays{j})=[interp_raw.(groups{i}).(assays{j})(:,1), fof];

       fbase=mean(data.(groups{i})(j).f_interp(30:250));
       fdelta=data.(groups{i})(j).f_interp-fbase;
       fof=fdelta/fbase;
       data.(groups{i})(j).dff_interp=fof;

       
   end

end

% interpolate tracking data
for j=1:numAssays
    XA=(1:maxTrack)/(maxTrack/210);
    B=tracking.(groups{1}).(assays{j})(:,4);
    XB=tracking.(groups{1}).(assays{j})(:,1);
    new_B=interp1(XB,B,XA);
    track_long_interp1.(assays{j})=[XA; new_B];

end

%% plot Ca transients through time
figure()
visOn=nan(); % set to time (in seconds) visual stimulus came on or nan() if no visual stimulus present
acvOn=90; % set to time (seconds) odor came on


%define visual stimulus start
x=dff.(groups{1}).(assays{1})(:,1);
if isnan(visOn)
    visBegin=length(dff.(groups{1}).(assays{1})(:,1));
else
    visBegin=min(find(x>visOn))-1;
end

for i = 1:numGroups
    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{i})));
   assays = fieldnames(dff.(groups{i}));

   dff_all=nan(size(dff.(groups{i}).(assays{j})(:,1),1),numAssays);
 
    subplot(numGroups+1,1,i)
    track_all=dff_all;
   for j = 1:numAssays

       subplot(numGroups+1,1,i)
       
       x=dff.(groups{i}).(assays{j})(:,1);fof= dff.(groups{i}).(assays{j})(:,2);
       

       plotter(x,fof,odorOn=acvOn,odorLength=60,LineWidth=0.2,alpha=0.2)
       hold on
       dff_all(:,j)=fof;
       track_all(:,j)=trackInterp.(groups{i}).(assays{j})(:,2);
       track_all(1:visBegin,j)=0;

        yline(0,'--','Color',[0.7 0.7 0.7])
       xline(visOn)
       xlim([0 210])
       ylabel('dF/F')
       box off
        set(gca, 'Color', 'none');
   end
   name=split(groups{i},'_');
   yline(0,'--','Color',[0.7 0.7 0.7])
   xline(visOn)
   if numAssays>1
        plotter(x,nanmean(dff_all'),odorOn=acvOn,odorLength=60,LineWidth=2)
   end
   legend(name{4})
   xlim([0 210])
   ylabel('dF/F')
   box off
set(gca, 'Color', 'none');


  
end

subplot(numGroups+1,1,i+1)
for k=1:j
    %figure(k)
    subplot(numGroups+1,1,i+1)

    plotter(x,track_all(:,k),odorOn=acvOn,odorLength=60,LineWidth=0.2,alpha=0.2)
    xlabel('time (sec)')
    xlim([0 210])
    ylabel('TI')
    xline(visOn)
    yline(0,'--','Color',[0.7 0.7 0.7])
    ylim([-1 1])
    box off
    set(gca, 'Color', 'none');
end

if numAssays>1
    
        plotter(x,nanmean(track_all'),odorOn=acvOn,odorLength=60,LineWidth=2)
end
xlabel('time (sec)')
xlim([0 210])
ylabel('TI')
xline(visOn)
yline(0,'--','Color',[0.7 0.7 0.7])
ylim([-1 1])
box off
set(gca, 'Color', 'none');
%% quantify df/f before, during, and after odor stimulus presentation

%select sectioning dependent on experimental design

%sections=[30 90;90 150; 150 210]; whichsection=2; %use for assays without visual stimulus
sections=[0 30;30 90; 90 150; 150 210]; whichsection=3; %use for assays with visual stimulus

figure()
for i = 1:numGroups
    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{i})));
   assays = fieldnames(dff.(groups{i}));

 
   assayMeans(i).means=nan(numAssays,size(sections,1));
   assayMeans(i).peaks=nan(numAssays,size(sections,1));

   subplot(1,numGroups,i)
   for j = 1:numAssays
        
        sectionmeans=nan(size(sections,1),1);
        for k=1:size(sections,1)
            from=sections(k,1)+5;
            to=sections(k,2)-5;

            range=find(dff.(groups{i}).(assays{j})(:,1)>from & dff.(groups{i}).(assays{j})(:,1)<=to);

            sectionmeans(k)=nanmean(dff.(groups{i}).(assays{j})(range,2));
            sectionpeaks(k)=max(dff.(groups{i}).(assays{j})(range,2));

        end
        assayMeans(i).means(j,:)=sectionmeans;
        assayMeans(i).peaks(j,:)=sectionpeaks;
   end
   x=repmat(1:size(sections,1),numAssays,1);

    scatter(x,assayMeans(i).means,'filled')
    hold on
    x2=repmat(1:size(sections,1),2,1)';
    x2(:,1)=x2(:,1)-0.25;
    x2(:,2)=x2(:,2)+0.25;
    y2=repmat(mean(assayMeans(i).means,1),2,1)';
    plot(x2',y2','k','LineWidth',2)
    plot(repmat(1:size(sections,1),numAssays,1)',assayMeans(i).means','Color',[0.7 0.7 0.7])
    xlim([0.5 size(sections,1)+0.5])
    ylim([-0.15 1.1])
    xticks(1:size(sections,1))
    %xticklabels({'pre','ACV','post','vis'}) %use for assays without visual stimulus
    xticklabels({'pre','vis','ACV+vis','vis'}) %use for assays with visual stimulus

    name=split(groups{i},'_');
    legend(name{4})
    ylabel('dF/F')
end




%% calculate and plot probability of males tracking in different stimulus conditions
%sections=[60 90;90 120]; %use for assays without visual stimulus
sections=[30 90; 90 150; 150 210]; %use for assays with visual stimulus
figure()

    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{i})));
   assays = fieldnames(dff.(groups{i}));

 
   assayMeans=nan(numAssays,size(sections,1));
   

       
       cutoff=0.05;
       for j = 1:numAssays
             %bind tracking index between 0 and 1
             ti=track_long_interp1.(assays{j})(2,:)';
             ti(ti<cutoff)=0;

             sectionmeans=nan(size(sections,1),1);
             for k=1:size(sections,1)
                 from=sections(k,1)+5;
                 to=sections(k,2)-5;

                  range=find(track_long_interp1.(assays{j})(1,:)>from & track_long_interp1.(assays{j})(1,:)<=to);
            
                  sectionmeans(k)=sum(ti(range))/length(range);
             end
             assayMeans(j,:)=sectionmeans;

       end
        
 
  
   x=repmat(1:size(sections,1),numAssays,1);

    scatter(x,assayMeans,'filled')
    hold on
    x2=repmat(1:size(sections,1),2,1)';
    x2(:,1)=x2(:,1)-0.25;
    x2(:,2)=x2(:,2)+0.25;
    y2=repmat(mean(assayMeans,1),2,1)';
    plot(x2',y2','k','LineWidth',2)
    plot(repmat(1:size(sections,1),numAssays,1)',assayMeans','Color',[0.7 0.7 0.7])
    xlim([0.5 size(sections,1)+0.5])
    ylim([0 1])
    xticks(1:size(sections,1))
    %xticklabels({'pre','ACV','ACV+vis','vis'}) %use for assays without visual stimulus
    xticklabels({'vis','ACV+vis','vis'}) %use for assays with visual stimulus
    ylabel('p(tracking)')


%% plot df/f against tracking index

%define ranges to plot dependent on experimental design
range=[1 210; 1 30; 31 90; 91 150; 151 210]; 
cr={'k','k','k','r','k','k'};
    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{i})));
   assays = fieldnames(dff.(groups{i}));
    
   for j=1:numAssays

       figure()
       p=1;
       for r=1:size(range,1)

           subplot(1,size(range,1),p)
           from=range(r,1);
           to=range(r,2);
           focus=find(dff.(groups{i}).(assays{j})(:,1)>from & dff.(groups{i}).(assays{j})(:,1)<to);
    
           fl=dff.(groups{i}).(assays{j})(focus,2);
           %option to normalize
           %fl=fl/max(fl);
           ti=trackInterp.(groups{i}).(assays{j})(focus,2);
           scatter(ti,fl,cr{r},'filled','MarkerFaceAlpha',.1)
           hold on

       
           xlim([-1 1])
           ylim([-0.3 3.5])
           [R,P]=corrcoef([ti fl]);
           text(-0.8,0.8, num2str(R(1,2)));
           text(-0.8,0.6,num2str(P(1,2)));
    
           % Fit a straight line (y = a*x + b)
            poly = polyfit(ti, fl, 1);   % 1 means linear fit
            yfit = polyval(poly, ti);
            plot(ti, yfit, 'r-', 'LineWidth', 2);
    
            if p==3
                xlabel('Tracking Index')
            end
            if p==1
                ylabel('dF/F')
            end


            titles={'all','pre','vis','vis+acv','vis'};

            title(titles{p})
            p=p+1;
          
            
       end
       
   end


%% plot dF/F recorded in LPC vs AOTu for fru>GCaMP and dsx>GCaMP assays


figure()

range=[1 210; 1 30; 31 90; 91 150; 151 210]; %use for assays with visual stimulus
cr={'k','k','k','r','k'};

%range=[1 210; 1 90; 91 150; 151 210]; %use for assays without visual stimulus
%cr={'k','k','r','k'};

p=1;    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{1})));
   assays = fieldnames(dff.(groups{1}));
    
   
   for r=1:size(range,1)
       subplot(numGroups,size(range,1),p)
       for j = 1:numAssays
            from=range(r,1);
            to=range(r,2);
           focus=find(dff.(groups{1}).(assays{j})(:,1)>from & dff.(groups{1}).(assays{j})(:,1)<to);
    
           fl_LPC=dff.(groups{1}).(assays{j})(focus,2);
           fl_AOTu=dff.(groups{2}).(assays{j})(focus,2);
           
           
           scatter(fl_AOTu,fl_LPC,cr{r},'filled','MarkerFaceAlpha',.1)
           hold on
       end
       plot([-10 10],[-10 10],'k')
       xlim([-0.3 0.4])
       ylim([-0.3 0.4])
       [R,P]=corrcoef([fl_AOTu fl_LPC]);
       text(-0.3,0.3, num2str(R(1,2)));
       text(-0.3,0.2,num2str(P(1,2)));
       p=p+1;
   end


%% calculate and plot p(of male tracking) over time
visOn=30; %define time visual stimulus comes on
acvOn=90; %define time odor comes on

x=track_long_interp1.(assays{1})(1,:);
visBegin=min(find(x>visOn))-1;


figure()
ti=nan(numAssays,size(track_long_interp1.(assays{j})(2,:),2));
cutoff=0.05;
for j = 1:numAssays
         ti(j,:)=track_long_interp1.(assays{j})(2,:)';
end
ti(ti<cutoff)=0;
ti(ti>=cutoff)=1;
tiTotal=sum(ti);
pTrack=tiTotal/numAssays;
pTrack(1:visBegin)=0;
pTrackSm=smoothdata(pTrack,"movmean",600);
plotter(x,pTrackSm,odorOn=acvOn,odorLength=60,LineWidth=1)
hold on
xline(visOn)
xlim([0 210])
xlabel("time (sec)")
ylabel("p(tracking)")

%% plot LC10 df/f vs positin of the visual stimulus


figure()

range=[1 210; 1 30; 31 90; 91 150; 151 210]; %vis30acv90

cr={'k','k','k','r','k','k'};
colM=[1 0 0; 0 0 0];
p=1;
for i = 1:numGroups
    
   %totals first 
   numAssays = length(fieldnames(dff.(groups{i})));
   assays = fieldnames(dff.(groups{i}));
    
   
   for r=1:size(range,1)
       subplot(numGroups,size(range,1),p)
       for j = 3 %:numAssays
           minAng=min(trackInterp.(groups{i}).(assays{j})(:,3));
           minAngAdd=minAng*-2;
           dffnorm=dff.(groups{i}).(assays{j})(:,2);
           dffnorm=dffnorm/max(dffnorm);
            visStimDir=zeros(1,size(dff.(groups{i}).(assays{j}),1));

            tmp=trackInterp.(groups{i}).(assays{j})(:,3);
            for l=2:size(trackInterp.(groups{i}).(assays{j})(:,3),1)
                if tmp(l)>tmp(l-1)
                    visStimDir(l)=1; %increasing = 1 decreasing =0
                end
                visStimDir(1)=visStimDir(2); %set first to second
            end

            from=range(r,1);
            to=range(r,2);
           focus=find(dff.(groups{i}).(assays{j})(:,1)>from & dff.(groups{i}).(assays{j})(:,1)<to);
    
           fl=dffnorm(focus);
           colRange=visStimDir(focus);
           colRange=colRange'+1;
           
           %option to normalize
           %fl=fl/max(fl);
           stim=trackInterp.(groups{i}).(assays{j})(focus,3);
           what=find(colRange==2);
           stim(what)=stim(what)*-1;
           stim(what)=stim(what)+minAngAdd;
           scatter(stim,fl,[],colM(colRange,:),'filled','MarkerFaceAlpha',.1)
           hold on

           saver(r).res=[stim fl];
       end
       xlim([-1.5 4])
       ylim([-0.3 1])

       p=p+1;
   end
end 
 
%% plot one individual assay

select=2; %select assay # to plot
ll=[0 210];
figure()
lastFrame=data.(groups{1})(select).track_size_raw;
x= (1:lastFrame)/(lastFrame/210);

subplot(4,1,1)
y=data.(groups{1})(select).track_femaleAngle_raw;
y=rad2deg(y);
y(1:find(x>=30,1))=nan();
plot(x,y,'k','LineWidth',1)
xlim(ll)
ylim([-100 100])
ylabel('Target position (deg)')
saver.target=y';
saver.time=x;

% Access current axes
ax = gca; 
% Remove all axes except the left y-axis
ax.Box = 'off';               % Turn off the surrounding box
ax.XColor = 'none';           % Hide the x-axis
ax.YColor = 'k';              % Keep the left y-axis visible
ax.Color='none';

subplot(4,1,2)
y=data.(groups{1})(select).track_hDiff_raw;
y=y*(numel(y)/210);
plotter(x,y,odorOn=90,odorLength=60,LineWidth=1)
xlim(ll)
ylim([-25 25])
ylabel('∆heading (rad*s-1)','FontName','Arial')

saver.heading=y';

% Access current axes
ax = gca; 
% Remove all axes except the left y-axis
ax.Box = 'off';               % Turn off the surrounding box
ax.XColor = 'none';           % Hide the x-axis
ax.YColor = 'k';              % Keep the left y-axis visible
ax.Color='none';

subplot(4,1,3)
y=data.(groups{1})(select).track_TI_raw;
y(1:find(x>=28,1))=nan();
yline(0, '--','Color',[.5 .5 .5])
hold on
plotter(x+2,y,odorOn=90,odorLength=60,LineWidth=1)
xlim(ll)
ylim([-1 1])
ylabel('Tracking index')

% Access current axes
ax = gca; 
% Remove all axes except the left y-axis
ax.Box = 'off';               % Turn off the surrounding box
ax.XColor = 'none';           % Hide the x-axis
ax.YColor = 'k';              % Keep the left y-axis visible
ax.Color='none';

subplot(4,1,4)
y=data.(groups{1})(select).track_maleAngVel_raw;
plotter(x,y,odorOn=90,odorLength=60,LineWidth=1)
xlim(ll)
ylabel('Angular velocity (rad*s-1)')

% Access current axes
ax = gca; 
% Remove all axes except the left y-axis
ax.Box = 'off';               % Turn off the surrounding box
ax.XColor = 'k';           % Hide the x-axis
ax.YColor = 'k';              % Keep the left y-axis visible
ax.Color='none';


%% functions

function [] = plotter(x,y,options)
    arguments
        x (1,:)
        y (1,:)
        options.odorOn (1,1) {mustBeNumeric} = 90
        options.odorLength (1,1) {mustBeNumeric} = 60
        options.LineWidth (1,1) {mustBeNumeric} = 0.5
        options.alpha (1,1) {mustBeNumeric} = 1
    end

    colsRGB=[0 0.6 1 options.alpha; 1 0.6 0 options.alpha];
    odorOff=options.odorOn+options.odorLength;
    start=0;
    stop=0;
    for i =2:length(x)
    
        if (start==0 && x(i)>options.odorOn)
            start=i;
    
        end
    
        if stop==0 && x(i)>odorOff
            stop=i;
        end
    end
    if isnan(options.odorOn)
        plot(x(1:end),y(1:end),'Color',[0 0 0 options.alpha],'LineWidth',options.LineWidth)
        hold on
    else
        plot(x(1:start),y(1:start),'Color',[0 0 0 options.alpha],'LineWidth',options.LineWidth)
        hold on
        plot(x(stop:end),y(stop:end),'Color',[0 0 0 options.alpha],'LineWidth',options.LineWidth)
        plot(x(start:stop),y(start:stop),'Color',colsRGB(2,:),'LineWidth', options.LineWidth)
    end
end





