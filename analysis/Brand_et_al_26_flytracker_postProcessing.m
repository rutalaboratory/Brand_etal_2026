% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% loads flytracker output, calculates basic parameters, identifies tracking
% and mating bouts and saves various parameters per fly in .mat file
% 
% to run select main folder of flytracker output for a single video
% requires a file called sex_index.txt to run. The file needs to list the
% sex of all flies based on their IDs in flytracker as a comma separated
% list, where 1 indicates a male and 0 indicates a female: 0,1,1,0,1,0
%
%
% the script is meant to be used interactively. How close male and female are during copulation
% differs between species - these parameters need to be individually identified
% by varying max_cop_mean and max_cop_sd


clear %clears workspace

% PARAMETERS
%
% initial parameter settings that help set distance thresholds to call interactions
% between flies. These vary based on resolution of individual videos (pixel per fly). 
% 
double=1; %set to 1 if 2 chambers were recorded at the same time (lower resolution than single assay)
quadrouple=0; %set to 1 if 4 chambers were recorded at the same time (low resolution per assay)

%set to one if you want to detect copulations 
score_mating=1;

%set to one if you want the script to plot parameter results (useful for trouble shooting)
print_int_o_time=1;
save_data=1;
minIntFrames=10; %number of frames used to solve conflicts in following etc. behavior %10 usually 30 for projector assay
gapClosingFrames=10; %max number of frames inbetween blocks to merge them %10 usually %30 for projector assay


max_interact_angle = pi/4; %the maximum angle subtended by the orientation of the interactor and the line connecting the two flies' center of area to count an interaction
max_follow_ang = deg2rad(15);
max_movement_angle = deg2rad(60);
max_dist = 4*avgLens; %4 for normal assays. 100 for projector assays
max_interact_dist = 1.5*avgLens;
min_vel = 0.5; %minimum velocity to score approaching and following behavior -> has to be empirically identified based on setup used
max_cop_mean = 60; %maximum mean interfly distance to ID mating
max_cop_sd = 10; %max sd interfly distance to ID mating


tic

fps = 30; %fps of tracked video

% LOAD DATA

nan = NaN();

% if there are multiple parts to a video, the script loads all matrices and concatenates them in order to analyze them together.
% Define a starting folder.
start_path = fullfile('/PATH/TO/FLYTRACKER/FOLDER');

% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ':');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

%initiate structure wherein flytracker output data is concatenated
flytracker = struct('x',[],'y',[],'length',[],'ori',[],'vel',[],'angvel',[],'minWingAng',[],'maxWingAng',[],'ar',[]);

%load necessary data from flytracker output
for k = 2:numberOfFolders
    
    thisFolder = listOfFolderNames{k}; %this is the correct way to read from a cell array
	fprintf('Processing folder %s\n', thisFolder);
    
    % Get trk files.
	filePattern = sprintf('%s/*track.mat', thisFolder);
	baseFileNames = dir(filePattern);
    baseFileNames = baseFileNames(~startsWith({baseFileNames.name}, '.'));
    numberOfMatFiles = length(baseFileNames);
    
    if numberOfMatFiles >= 1
	% Go through all those mat files (probably only one per folder but who knows).
        for f = 1 : numberOfMatFiles
            fullFileName = fullfile(thisFolder, baseFileNames(f).name);
            fprintf('     Processing trk file %s\n', fullFileName);

            load(fullFileName)
            %attach the data in the trk struct you need into one matrix
            flytracker.x = [flytracker.x; trk.data(:,:,1)'];
            flytracker.y = [flytracker.y; trk.data(:,:,2)'];
            flytracker.length = [flytracker.length; trk.data(:,:,4)'];
            flytracker.ori = [flytracker.ori; trk.data(:,:,3)'];

            featFileName = strrep(fullFileName,'track.mat','feat.mat');
            fprintf('     Processing feat file %s\n', featFileName);

            load(featFileName)
            flytracker.vel = [flytracker.vel; feat.data(:,:,1)'];
            flytracker.angvel = [flytracker.angvel; feat.data(:,:,2)'];
            flytracker.minWingAng = [flytracker.minWingAng; feat.data(:,:,3)'];
            flytracker.maxWingAng = [flytracker.maxWingAng; feat.data(:,:,4)'];
            flytracker.ar = [flytracker.ar; feat.data(:,:,6)'];

        end


	else
		fprintf('     Folder %s has no track files in it.\n', thisFolder);
	end
    
end
fprintf('Finished loading data\n\n');
%how many flies are in the assay?
numflies = size(flytracker.x,2); 

%load sex information
sex_on=0;

fprintf('looking for sex index\t');

sex_file = append(topLevelFolder,'/',"sex_index.txt");

if isfile(sex_file)
    fprintf('-- found it\n\n');

    sex = readmatrix(sex_file,'Delimiter',',');
    sex_on =1;
else
    fprintf("-- couldn't find it - make sure to provide a file called sex_index.txt\n\n");
end



%fix swaps in fly orientation (happens when flies walk backwards
%flytracker infers the flies front based on walking direction
ori_swap_ids={};
for num = 1:numflies
    
    mean_ar = nanmean(flytracker.ar(:,num));
    focfly = flytracker.ori(:,num);
    swp_tmps=0;
    
    last_potential_swap=nan();
    
    i=2;
    while (i<=size(flytracker.ori,1))
        swp=0;
        
        deltaOri=0;
        if ~isnan(focfly(i-1))
            deltaOri = angle_difference(focfly(i),focfly(i-1)); 
        else
            if i==2 && isnan(focfly(1))
                ii=1;
                while isnan(focfly(i+ii))
                   
                    ii=ii+1;
                    
                end
                i=i+ii;
            end
            z=2;
            while isnan(focfly(i-z))
               z=z+1; 
               
               if z == i-1
                  
                  deltaOri=0;
                  break
               end
            end
            
            if z<i-1 && i > 2
                deltaOri = angle_difference(focfly(i),focfly(i-z-1));
            else
                deltaOri=0;
            end
        end
       %fix angle if swapped
       if deltaOri>deg2rad(120) && deltaOri<deg2rad(340) %second term makes sure not to 'correct' when fly moves from + to - part of the Einheitskreis
           swp=1;
           
           min_index=i-3*fps;%make sure min_index is not out of bounds
           if min_index<1
               min_index=1;
           end
               
           if max(flytracker.maxWingAng(min_index:i,num))>1.5 || max(flytracker.vel(min_index:i,num))>40 || max(flytracker.angvel(min_index:i,num))>=20  || i-last_potential_swap<6*fps || min(flytracker.ar(min_index:i,num))<mean_ar*0.7
               %this term excludes flying individuals or individuals that are not faithfully detected (small axis ratio) - these often swap oris but then stabilize right after
               %also excludes swaps that happen within 3 seconds after a
               %potential swap that was not executed due to the first 3
               %terms.
                swp=0; %Original
               %nothing to be done here with exception of ongoing flip
                if swp_tmps(end)+1 == i
                    swp=1;

                else                    
                    if deltaOri>deg2rad(160) && max(flytracker.maxWingAng(min_index:i,num))>1.5 && max(flytracker.vel(min_index:i,num))<5 %this might happen if individuals are clustered but not moving much
                        swp=1;
                    else
                        last_potential_swap=i;
                    end
                    
                end
                

           end
           if swp
               
               tmp = pi-abs(focfly(i));
               if focfly(i)>0 %if the previous (and incorrect) angle was positive, the opposite has to be negative
                   tmp = tmp * -1; 

               end
               focfly(i) = tmp; 

               %delete erroneous wing detection during wrong orientation
               flytracker.minWingAng(i,num) = nan;
               flytracker.maxWingAng(i,num) = nan;
               swp_tmps=[swp_tmps i];
           end
       end
       i=i+1;
    end
    
      
    flytracker.ori(:,num) = focfly;
    ori_swap_ids{num,1}=swp_tmps;
end


%get mean length for each fly and for all flies
lens_means = zeros(numflies,1);
for ind = 1:numflies
    tmp_m = nanmean(flytracker.length(:,ind));
    lens_means(ind,1) = tmp_m;
end

avgLens = mean(lens_means); %this size will change from assay to assay dependent on how many flies of which sex are present, how many copulations etc. but does not change performance
 

numcomp = sum (numflies-1:-1:1);
sec = 11;
wsize=fps*sec;  %window size 
jump=wsize/2;    %window slide
seq = 1:jump:size(flytracker.x,1)-wsize +1;

res = zeros(size (seq,2), numcomp*2 ); %matrix to fill
lens= zeros(size (seq,2), numflies);


interaction_map = zeros(size (flytracker.x,1),numflies);
approach_map = interaction_map; %first row will always remain empty because it is based on change in position between frames.
mating_map = interaction_map;
wing_map = interaction_map;


%initiate cell arrays
interaction_cell = cell(numflies,1);
for c = 1:numflies
    interaction_cell{c,1}=interaction_map;    
end
approach_cell = interaction_cell;
rel_position_cell = interaction_cell;
rel_ori_cell = interaction_cell;
rel_walking_dir_cell = interaction_cell;
euc_dist_cell = interaction_cell;
mating_cell=interaction_cell;

walking_direction = zeros(size(flytracker.x,1)-1,numflies);
wd_check = 0;


if double
    max_cop_mean=max_cop_mean/2; 
    max_cop_sd = max_cop_sd/2;
end

if quadrouple
    max_cop_mean=max_cop_mean/3; 
    max_cop_sd = max_cop_sd/3;
end

it = 1; %increases by 2
it2 = 1; %increases by 1
    for fly1 = 1:numflies-1
        
       
        
        for fly2 = fly1+1:numflies %this ensures that each combination is only looked at once
            
            
            
            dir_diff = zeros(size(walking_direction,1),2); %saves difference in orientation angles to angle connecting 2 flies. helps identifying approaching and following behavior.

            
            %%%%this part of the code calculates if 
            %%%%1. focal fly is oriented towards and walks target fly
            
         
            % extract x and y coordinates of two focal flies
            f1 = [flytracker.x(:,fly1),flytracker.y(:,fly1)];
            f2 = [flytracker.x(:,fly2),flytracker.y(:,fly2)];

            % calculate euclidean distance between them for all frames
            fly12 = sqrt((f1(:, 1) - f2(:, 1)) .^ 2 + (f1(:, 2) - f2(:, 2)) .^ 2);
            euc_dist_cell{fly1,1}(:,fly2)=fly12;
            euc_dist_cell{fly2,1}(:,fly1)=fly12;
            
            % identify tracking
            % 1. identify if focal fly is pointing towards another fly 
                % 1.1 needs to be within certain distance
                % 1.2 flies heading needs to be within certain degree angle
                % from other fly

           
            %check if 2 flies are close enough to each other to potentially interact
                
            close_enough = fly12 <= max_dist;
            close_enough_interaction = fly12 <=max_interact_dist;
            %check if focal fly is oriented towards other fly. Then
            %switch focal fly and check the reverse
            
            %this calculates the angle of the line between the two flies on
            %the unit circle
            rad2ref1 = atan2( f2(:,2)-f1(:,2), f2(:,1)-f1(:,1) ) * -1;   % multiply by one bc y is mirrored in our setup     
            
          % translate angle to complement on unit circle (this does take into account +- rad)
            rad1ref2 = pi - abs(rad2ref1);
            neg = find(rad2ref1>0);
            rad1ref2(neg) = rad1ref2(neg) * -1 ; 
            
            %make a copy of the arrays to save position angles. have to do
            %it, because rad2ref1 and rad1ref2 are going to be modified to
            %represent the orientation relative to position angle.
            pos_rad21 = rad2ref1;
            pos_rad12 = rad1ref2;
            
            %this is the actual angle we are looking for. The one that
            %describes the angle between the angle of the two flies and the
            %orientation of the focal fly
            %this means, we want the DIFFERENCE of the angles IF both have
            %the same sign
            %but the SUM of the absolute values of the angles IF they have
            %different signs.
            for x = 1:size(rad2ref1,1)
                
                if x>1
                    if wd_check < fly1 %makes sure to only calculate these once per fly
                       walking_direction(x-1,fly1) = atan2( f1(x,2)-f1(x-1,2), f1(x,1)-f1(x-1,1) ) * -1; %b/c y is mirrored

                    end
                    if wd_check < fly2
                       walking_direction(x-1,fly2) = atan2( f2(x,2)-f2(x-1,2), f2(x,1)-f2(x-1,1) ) * -1; %same

                    end
                end

                
                if close_enough(x) %%this part deletes a lot of angles calculated. 
                
                    
                    rad2ref1(x) = angle_difference(rad2ref1(x),flytracker.ori(x,fly1));
                    rad1ref2(x) = angle_difference(rad1ref2(x),flytracker.ori(x,fly2));
                    
                     
                %%%calculate the direction of movement of each fly. This
                %%%helps me extract following behavior (I'm doing it in
                %%%this for loop so I don't have to loop twice)

                if x>1
                    
                    players = [fly1; fly2];
                    for ind = 1:2
                        if ind ==1
                            target = fly2;
                            position = pos_rad21(x);
                        else
                            position = pos_rad12(x);
                            target = fly1;
                        end
                      

                         dir_diff(x-1,ind) = angle_difference(walking_direction(x-1,players(ind)),position);

                        %is fly1 walking towards fly2 in an acceptable
                        %angle between walking dir and line connecting the two flies and vice versa?
                        
                        if close_enough_interaction(x) %this adjusts max angle to distance between flies in a two-step system
                            acute_angle = max_interact_angle;
                        else
                            acute_angle = max_follow_ang;
                        end

                        if dir_diff(x-1,ind)<=acute_angle && flytracker.vel(x,players(ind)) >= min_vel
                            
                            if angle_difference(walking_direction(x-1,players(ind)),flytracker.ori(x,players(ind)))<=max_movement_angle %this checks if the fly is walking forwards by requiring that the difference between the orientation angle and the walking direction is within a certain angle 
                                
                                %approach_map(x,players(ind)) = target;
                                approach_cell{players(ind),1}(x,target)=target;
                                
                            end
                        end
                        
                        
                        
                    end
                end
                
                
                    
                    
                    
                    
                else
                    rad2ref1(x) = nan;
                    rad1ref2(x) = nan;
                end
            
                %save interaction if individuals are close enough and
                %interactor is at an angle < max_angle
                
                if rad2ref1(x)<max_interact_angle && close_enough_interaction(x)
                    %interaction_map(x,fly1) = fly2;
                    interaction_cell{fly1,1}(x,fly2)=fly2;

                end
                if rad1ref2(x)<max_interact_angle && close_enough_interaction(x)
                    %interaction_map(x,fly2) = fly1;
                    interaction_cell{fly2,1}(x,fly1)=fly1;
                end

            
            end
            
            if wd_check<fly2
                wd_check = fly2;
            end
            
            %save orientations and positions for subsequent map
            %constructions and conflict solving.
            rel_position_cell{fly1,1}(:,fly2)=pos_rad21;
            rel_position_cell{fly2,1}(:,fly1)=pos_rad12;

            rel_ori_cell{fly1,1}(:,fly2)=rad2ref1;
            rel_ori_cell{fly2,1}(:,fly1)=rad1ref2;

            rel_walking_dir_cell{fly1,1}(2:size(flytracker.x,1),fly2)=dir_diff(:,1);
            rel_walking_dir_cell{fly2,1}(2:size(flytracker.x,1),fly1)=dir_diff(:,2);

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%sliding window mean distance to identify mating pairs%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if score_mating
                c = 1;
                for frame = seq
                win = frame:(frame + wsize - 1);
                pdistwin = fly12(win); %current window
                res(c,it) = nanmean(pdistwin); %mean interfly distance in window
                res(c,it+1) = nanstd(pdistwin); %sd of interfly distance in window         

                
               hetero=0;
               if sex_on 
                    if sex(fly1) ~= sex(fly2)
                        hetero=1;
                    end
                else
                    if abs(lens_means(fly1)-lens_means(fly2))>std(lens_means)
                       hetero=1; 
                    end
                end


                if hetero && res(c,it)<max_cop_mean && res(c,it+1)<max_cop_sd 


                    mating_cell{fly1,1}(win,fly2)=fly2;
                    mating_cell{fly2,1}(win,fly1)=fly1;

                end

                    c = c+1;
                end

            end
          xmax = size(res,1);

          it = it+2;
          it2 = it2+1;
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%CORRECT CELLS WHEN THERE ARE CONFLICTS%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    mating_map=conflict_resolver(mating_cell,fps*60*5,euc_dist_cell); %5 mins - change if copulations are not correctly determined
    mating_map_bup=mating_map;
[mating_blocks,mating_map] = blocker(mating_map,fps*60*5);

interaction_map=conflict_resolver(interaction_cell,minIntFrames,euc_dist_cell); 
approach_map=conflict_resolver(approach_cell,minIntFrames,euc_dist_cell);


toc

mating_blocks_bup=mating_blocks;
for i =1:size(mating_blocks,1)
    %fprintf( 'working on animal %i\n', i)
   deleter = [];
    if size(mating_blocks{i,1})%if potential matings were recorded for that individual
        for j = 1:size(mating_blocks{i,1},1)
            
           
            
            to = mating_blocks{i,1}(j,1);
            from = to-fps; %one second
            if from <=0
                from=1;
            end
            interactee = mating_blocks{i,1}(j,3);
            duration=mating_blocks{i,1}(j,4);
                %fprintf('\tOn mating block %i with interactee %i\n',j,interactee)
             
                %fprintf('\t\tWindow from %i to %i\n',from,to)
            
            hetero=0;
            if sex_on 
                if sex(i) ~= sex(interactee)
                    hetero=1;
                end
            else
                if abs(lens_means(i)-lens_means(interactee))>std(lens_means) %this makes sure only mating between male and female is going to be looked at for mating
                   hetero=1; 
                end
            end
            
            
            % this part helps to remove wrongly detected copulations - this
            % happens sometimes if the flies stand close to one another. If
            % potential cases are identified the script puts out a mating
            % detection warning. use the `deleter` array to remove wrong
            % copulations.

            if duration>fps*60*5 %vary if mating detection performance is bad
                if hetero
                    if sum(interaction_cell{i,1}(from:to,interactee)) || sum(interaction_cell{interactee,1}(from:to,i))
                        %fprintf('\t\t\tInteraction detected\n')
                        %then nothing needs to be deleted
                        if nanmean(rel_ori_cell{i,1}(from:to,interactee))>2.8 || nanmean(rel_ori_cell{i,1}(from:to,interactee))<0.3 %this checks if individuals were oriented front to back in the past 10 frames prior to copulation
                            if nanmean(rel_ori_cell{i,1}(from:to,interactee))>2.8 && nanmean(rel_ori_cell{interactee,1}(from:to,i))<0.3
                           %     fprintf('\t\t\tfly%i is interactee, fly%i is interactor\n',i,interactee)
                            elseif  nanmean(rel_ori_cell{i,1}(from:to,interactee))<0.3 && nanmean(rel_ori_cell{interactee,1}(from:to,i))>2.8
                              %  fprintf('\t\t\tfly%i is interactee, fly%i is interactor\n',interactee,i)

                            else
                               %deleter = [deleter j];

                               fprintf('\t\t\tPotential mating detection error: not head to butt 1 - discard\n')

                            end
                        else

                            if nanmean(flytracker.vel(mating_blocks{i,1}(j,1):mating_blocks{i,1}(j,2)))>1.5
                                   fprintf('\t\t\tPotential mating detection error: special case 1.1 in mating map - check it!\n')
                                    %deleter = [deleter j];
                            else
                                fprintf('\t\t\tPotential mating detection error: special case 1.2 in mating map - check it!\n')
                                    %deleter = [deleter j];
                            end
                           fprintf('\t\t\tPotential mating detection error: not head to butt 2 - discard\n')
                           %deleter = [deleter j];
                        end
                    else
                        %this means there is no recorded interaction within 10
                        %frames before potential mating block starts

                        if nanmean(flytracker.vel(mating_blocks{i,1}(j,1):mating_blocks{i,1}(j,2)))>1.5 && length(flytracker.vel(mating_blocks{i,1}(j,1):mating_blocks{i,1}(j,2)))>fps*60*15
                            fprintf('\t\t\tPotential mating detection error: special case 2.1 in mating map - check it!\n')
                            
                        else
                            fprintf('\t\t\tPotential mating detection error: special case 2.1 in mating map - check it!\n')
                                %deleter = [deleter j];
                        end
                    end
                else
                    deleter = [deleter j];
                end
            else
                deleter = [deleter j];
            end
        end
        
    end
    
     if size(deleter,1)
        
         for del =length(deleter):-1:1 %delete from last to first
             %fprintf('\t\t\t\tWill delete mating block number %i\n',deleter(del))
             %fprintf('from %i to %i set to zero\n',mating_blocks{i,1}(deleter(del),1),mating_blocks{i,1}(deleter(del),2))
            
              %translate maps back - set from NaN to 0
            rebound_start=mating_blocks{i,1}(deleter(del),1);
            rebound_end=mating_blocks{i,1}(deleter(del),2);
            rebound_interactee=mating_blocks{i,1}(deleter(del),3);
            
            interaction_map(rebound_start:rebound_end,i)= interaction_cell{i,1}(rebound_start:rebound_end,rebound_interactee) ;
            
            approach_map(rebound_start:rebound_end,i)= approach_cell{i,1}(rebound_start:rebound_end,rebound_interactee) ;
            
            mating_map(rebound_start:rebound_end,i)= 0;
            
             mating_blocks{i,1}(deleter(del),:) =[];
             
         end
    end
    
end


%remove all wing parameters from mating individuals, bc they are msotly
%rubbish (flytracker is not great at detecting those)
%also remove interactions and approach calls, because mating individuals
%are defined as not approaching or interacting with flies other than the
%mating partner
%However, incoming interactions are still valid.

for line = 1: size(mating_map,1)
    for indiv = 1:numflies
       if mating_map(line,indiv)
          
          flytracker.minWingAng(line,indiv) = nan;
          flytracker.maxWingAng(line,indiv) = nan;
          interaction_map(line,indiv) = nan;
          approach_map(line,indiv) = nan;
       end
    end
end


%analyze wings
 lag = 10; %10
 threshold = 1.5; %1.5
 influence = 0;
 
for foc_fly = 1:numflies  
    
    wing_map(:,foc_fly) = detect_peaks(flytracker.maxWingAng(:,foc_fly),lag,threshold,influence);

    for i = 1:size(wing_map(:,1),1)
       if wing_map(i,foc_fly)>0
          if approach_map(i,foc_fly)>0
              wing_map(i,foc_fly)=approach_map(i,foc_fly);
          elseif interaction_map(i,foc_fly)>0
              wing_map(i,foc_fly)=interaction_map(i,foc_fly);
          else
              wing_map(i,foc_fly)=0;
          end
       end
    end
    
end


%closing gaps
pre_interaction_map2=interaction_map;
[interaction_blocks,interaction_map] = blocker(interaction_map,gapClosingFrames);
[approach_blocks,approach_map] = blocker(approach_map,gapClosingFrames);
[wing_blocks,~] = blocker(wing_map,0);

%extract courtship - requires sex information
%remove called interactions if focal fly is immobile

courtship_map=zeros(size(interaction_map,1),numflies);
female_following_map=courtship_map;

%singing detection - coded it but ended up not using it because wing
%detection is terrible in flytracker
prepost_max = fps; %seconds pre and post interaction to look for singing
min_int_vel=0.5; 
min_int_angvel=0.5; 

for myfly = 1:numflies
    
    focsex=sex(myfly);
    
   kill_blocks=[];
    for i = 1:size(interaction_blocks{myfly,1},1)
       
        now=interaction_blocks{myfly,1}(i,:); 
        
        tofly=now(3);
        tosex=sex(tofly);
        vel_m=nanmean(flytracker.vel(now(1):now(2),myfly));
        angvel_m=nanmean(flytracker.angvel(now(1):now(2),myfly));
        
        call=0;
        if vel_m>min_int_vel || angvel_m>min_int_angvel
            call=1;
        else
            prepost=round(now(4)/10);
            
            if prepost>prepost_max
               prepost=prepost_max; 
            end
            
            
            from=now(1)-prepost;
            to=now(2)+prepost;
            if from<1
               from=1; 
            end
            if to>size(wing_map,1)
               to=size(wing_map,1); 
            end
            if sum(wing_map(from:to,myfly))
                
                if sum(wing_map(from:now(1)+prepost,myfly)) && sum(wing_map(now(2)-prepost:to,myfly))
                    %does he sing before and after mostly immobile
                    %'interaction'?
                    call=1;
                    
                else
                    %special case -> male lost attention during sitting.
                    %Call half of the bout courtship up to max of mean
                    %courtship bout length of species to be conservative. 
                    call=2;
                end
                
            end
        end
        
        if call
            
            range = now(1):now(2);
            if call ==2

                maxi=fps*30; 
                
                if now(4)/2 > maxi
                   range = now(1):now(1)+maxi-1;
                else
                    range = now(1):now(1)+(now(4)/2)-1;
                end
                
                %correct interaction_block accordingly
                interaction_map(max(range)+1:now(2),myfly)=0;
                interaction_blocks{myfly,1}(i,2)=max(range); 
                interaction_blocks{myfly,1}(i,4)=(max(range)-now(1))+1;
                

                
            end
            
            if focsex==1 && tosex==0

                courtship_map(range,myfly)=tofly;
            else
                if focsex ==1 && tosex == 1
                   %%check if male is tracking other male. 
                   
                   if nnz(approach_map(range,myfly)>0)>=now(4)/4
                       %%if approach is called for 1/4 of interaction call
                       %%it tracking (can still be altered later)
                       
                       %later put in restriction to myfly targeting toflies
                       %back
                       
                       courtship_map(range,myfly)=tofly;
                   else
                       %look for wing extensions
                       
                       if sum(wing_map(range,myfly))
                           %call courtship for now - build in a subroutine
                           %that can distinguish aggression from singing
                           courtship_map(range,myfly)=tofly;
                       end
                   
                   end
                   
                else
                    if focsex == 0
                        if nnz(approach_map(range,myfly)>0)>=now(4)/2
                            
                            female_following_map(range,myfly)=tofly;
                            
                            %also add following outside of interaction to
                            %map;
                            
                        end
                    end
                    
                end

            end
        
            
        else
            %if call=0, which means the focal individual (myfly) is
            %motionless for too long (stopped interacting close to other
            %fly)
            %delete interaction block

            kill_blocks=[kill_blocks i];
        end
        
    end
    
    %now delete deleted interaction blocks
    for i = length(kill_blocks):-1:1
        focal=interaction_blocks{myfly,1}(kill_blocks(i),:);
        
        interaction_map(focal(1):focal(2),myfly)=0;
        
        interaction_blocks{myfly,1}(kill_blocks(i),:)=[];
    end
    
    
    
end

%plot results to check performance
if print_int_o_time
    num_frames = size(flytracker.x,1);
    block_printer(interaction_blocks,numflies,0,'b',num_frames,sex,10);
    block_printer(approach_blocks,numflies,0.1,'r',num_frames,sex,10);
    block_printer(wing_blocks,numflies,0.2,'k',num_frames,sex,1);
    block_printer(mating_blocks,numflies,0.3,'m',num_frames,sex,10);
end



toc



% save data output
if save_data
    results = struct('flies',[],'sex',[],'fly_lengths',[],'x',[],'y',[],...
    'approach_map',[],'interaction_map',[],'wing_map',[],...
    'mating_map',[],'courtship_map',[],'approach_blocks',[],'interaction_blocks',[],'mating_blocks',[],...
    'rel_ori',[],'rel_position',[],'rel_walking_direction',[],'vel',[],'dist',[]);

    results.flies=1:numflies;
    results.sex=sex;
    results.fly_lengths=lens_means;
    results.x=flytracker.x;
    results.y=flytracker.y;
    results.approach_map=approach_map;
    results.interaction_map=interaction_map;
    results.wing_map=wing_map;
    results.mating_map=mating_map;
    results.courtship_map=courtship_map;
    results.approach_blocks=approach_blocks;
    results.interaction_blocks=interaction_blocks;
    results.mating_blocks=mating_blocks;
    results.rel_ori=rel_ori_cell;
    results.rel_position=rel_position_cell;
    results.rel_walking_direction=rel_walking_dir_cell;
    results.vel=flytracker.vel;
    results.dist=euc_dist_cell;
    results.fps=fps;
    
    [filepath,name,ext] = fileparts(topLevelFolder);
    filename = strcat(topLevelFolder,'/',name,'.mat');
    save(filename,'results');
end


%ci
ii=zeros(1,numflies); %interaction index
ii_m = ii;
ii_f = ii;
ii_firstMating = ii;
ai=ii;
ai_m = ii;
ai_f=ii;

ai_map=zeros(size(approach_map,1),numflies);
focal_map=ai_map;

for i = 1:numflies
    from=1;
    to= size(approach_map,1);
    to2=to;
    

    
    %merge interaction and approach map to 'holistic interaction map'
    focal_map(:,i)=interaction_map(:,i);
    
    for j = 1: to
       
          if(approach_map(j,i)>0)
              if focal_map(j,i)==0
                focal_map(j,i)=approach_map(j,i); 
              end
             ai_map(j,approach_map(j,i))=i;
          end
          
          if(interaction_map(j,i)>0)
              ai_map(j,interaction_map(j,i))=i;
          end
       
    end
    
    %ii only up until first mating
    if ~isempty(mating_blocks{i,1})
        to2 = mating_blocks{i,1}(1,1)-1;
    end
    total_interactions = nnz(focal_map(from:to,i)>0);
    
    
    ii(i)=total_interactions/to;
    ii_firstMating(i)=nnz(focal_map(from:to2,i)>0)/to2;
    
    ii_m(i)=nnz(sex(focal_map(focal_map(from:to,i)>0,i))==1)/to;
    ii_f(i)=nnz(sex(focal_map(focal_map(from:to,i)>0,i))==0)/to;

end
for j = 1:numflies
    ai(j)=nnz(ai_map(:,j)>0)/size(approach_map,1);
    ai_m(j)= nnz(sex(ai_map(ai_map(:,j)>0,j))==1)/size(approach_map,1);
    ai_f(j)= nnz(sex(ai_map(ai_map(:,j)>0,j))==0)/size(approach_map,1);
end


%make fraction interacting map



%mean velocity per indiviual over time

win_size=fps*30; %thirty second windows

steps = 1:win_size:size(flytracker.vel,1);

vel_bins = zeros(length(steps)-1,numflies);

for i = 2:length(steps)
   
    
    for j = 1:numflies
              
        vel_bins(i-1,j) = nanmean( flytracker.vel( steps(i)-win_size:steps(i)-1, j));

    end
    
end

%interactions per indivual
int_per_ind=zeros(numflies,1);
int_to_male=int_per_ind;
int_to_female=int_per_ind;

for block = 1:numflies
    
int_per_ind(block) = nnz(interaction_blocks{block,1}(:,4)>9); %blocks of 10 consequtive frames or more
    
    if sex_on
        %interactions seperated by sex
        int_to_male(block)=nnz(sex(interaction_blocks{block,1}(interaction_blocks{block,1}(:,4)>9,3))==1);
        int_to_female(block)=nnz(sex(interaction_blocks{block,1}(interaction_blocks{block,1}(:,4)>9,3))==0);
    end
    

    
    
end

figure(numflies+2)

%set standard colors
if sex_on
    newcolors = cell(1,numflies);

    for i = 1:numflies
       if sex(i) == 1
           
           newcolors{1,i}='#FF9900';
       else
           newcolors{1,i}='#00AA00';
           
       end
    end
    
    colororder(newcolors)
end



subplot(4,1,1)
b1 = bar([ii_firstMating;ii;ii_f;ii_m;ai;ai_f;ai_m]');
b1(1).FaceColor = '#000000';
b1(2).FaceColor = '#DD0000'; %makes first bars red
b1(3).FaceColor = '#00AA00';
b1(4).FaceColor = '#FF9900';
b1(5).FaceColor = '#0033FF'; %makes second bars blue
b1(6).FaceColor = '#00AA00';
b1(7).FaceColor = '#FF9900';

ylim([0 1]);
ylabel('ii and ai per individual');
xlabel('flies');

subplot(4,1,2)

ii_ratio_m=zeros(length(ii_f),1);
ii_ratio_f=ii_ratio_m;
ii_ratio=ii_ratio_m;

ai_ratio_m=ii_ratio_m;
ai_ratio_f=ii_ratio_m;
ai_ratio=ai_ratio_m;

for i =1:length(ii_f)
    ii_ratio_f(i)=ii_f(i)/ii(i);
    ii_ratio_m(i)=ii_m(i)/ii(i);
    ii_ratio(i)=ii(i)/(ii(i)+ai(i));
    
    ai_ratio_f(i)=ai_f(i)/ai(i);
    ai_ratio_m(i)=ai_m(i)/ai(i);
    ai_ratio(i)=ai(i)/(ii(i)+ai(i));
end

b3=bar([ii_ratio,ii_ratio_f,ii_ratio_m,ai_ratio,ai_ratio_f,ai_ratio_m]);
b3(1).FaceColor = '#DD0000'; %makes first bars red
b3(2).FaceColor = '#00AA00';
b3(3).FaceColor = '#FF9900';
b3(4).FaceColor = '#0033FF'; %makes second bars blue
b3(5).FaceColor = '#00AA00';
b3(6).FaceColor = '#FF9900';
ylabel('ii and ai female/male ratio');
xlabel('flies');
ylim([0 1]);

subplot(4,1,3)

b2=bar([int_per_ind,int_to_female,int_to_male]);

b2(1).FaceColor = '#880088'; %all interactions
b2(2).FaceColor = '#00AA00'; %interactions targeting females
b2(3).FaceColor = '#FF9900'; %interactions targeting males
ylabel('number of interactions');
xlabel('flies');
ylim( [0 3000] );




subplot(4,1,4)
plot(vel_bins)
hold on
plot(mean(vel_bins,2),'Linewidth',2,'Color','k')
if sex_on
    plot(mean(vel_bins(:,sex==1),2),'Linewidth',2,'Color','#FF9900'); %males
    plot(mean(vel_bins(:,sex==0),2),'Linewidth',2,'Color','#00AA00'); %females
end
hold off

ylim([0 15]);

ylabel('velocity (mm/s)');
xlabel('time in 30 second bins');






win_size=fps; %one second windows

steps = 1:win_size:size(flytracker.vel,1);

m_m=zeros(length(steps)-1,1);
m_f=m_m;
f_m=m_m;
f_f=m_m;

for i = 2:length(steps)
    n_males=sum(sex==1);
    n_females=sum(sex==0);
    
    males_to_females=0;
    females_to_males=0;
    males_to_males=0;
    females_to_females=0;
    
    for j = 1:numflies
        to_male=0;
        to_female=0;
        nans=0;
                
        uni=unique(focal_map(steps(i)-win_size:steps(i)-1,j));
        
        for each = 1:length(uni)
            %fprintf('yo!\n')
            if isnan(uni(each))
                %fprintf('yo!\n')

                nans=nans+1;
                
            else
                if uni(each)>0
                    if sex(uni(each))==1
                       to_male =1;
                    else 
                        to_female =1;
                    end
                end
            end
            
        end
        %fprintf('num nans: %d\n',nans);
        if nans == win_size
           %fprintf('yo! ind:%d win: %d - %d\n',j,steps(i)-win_size,steps(i));
            if sex(j)==1
               n_males=n_males-1;
            else
                n_females=n_females-1;
            end
        end
        
        if sex(j)==1
           males_to_females=males_to_females+to_female;
           males_to_males=males_to_males+to_male;
        else
            females_to_females=females_to_females+to_female;
            females_to_males=females_to_males+to_male;
        end

    end
    
    m_m(i-1)=males_to_males/n_males;
    m_f(i-1)=males_to_females/n_males;
    f_m(i-1)=females_to_males/n_females;
    f_f(i-1)=females_to_females/n_females;
        
    
end

%sliding window over 1 sec bins
w_size=10; %window size
slide=5;
step=10; %end of current window;

mm_sw=[];
mf_sw=[];
fm_sw=[];
ff_sw=[];
i=1;
while step <length(m_m)
    
    mm_sw(i)=nanmean(m_m((step-w_size)+1:step));
    mf_sw(i)=nanmean(m_f((step-w_size)+1:step));
    fm_sw(i)=nanmean(f_m((step-w_size)+1:step));
    ff_sw(i)=nanmean(f_f((step-w_size)+1:step));
    
    step=step+slide;
    i=i+1;
    
end

figure(numflies+3)


subplot(5,2,[1,2])
x=1:length(mm_sw);
plot(x,ff_sw,'g')
ylim([-0.05 1.05])
ylabel('fraction of individuals');
xlabel('time in 10 second bins');

subplot(5,2,[3,4])
plot(x,fm_sw,'k')
ylim([-0.05 1.05])
ylabel('fraction of individuals');
xlabel('time in 10 second bins');

subplot(5,2,[5,6])
plot(x,mm_sw,'r')
ylim([-0.05 1.05])
ylabel('fraction of individuals');
xlabel('time in 10 second bins');

subplot(5,2,[7,8])
plot(x,mf_sw,'b')
ylim([-0.05 1.05])
ylabel('fraction of individuals');
xlabel('time in 10 second bins');


subplot(5,2,9)
ls1=mf_sw'\mm_sw';
yCalc1=ls1*mf_sw;
scatter(mf_sw,mm_sw,'filled')
hline=refline(1,0);
hline.Color ='k';
hline.LineStyle='--';
hold on
plot(mf_sw,yCalc1)
ylabel('male to male');
xlabel('male to female');
ylim([0 1])
xlim([0 1])
hold off
[rm,pm]=corr([mf_sw; mm_sw]');

subplot(5,2,10)
ls1=ff_sw'\fm_sw';
yCalc1=ls1*ff_sw;

scatter(ff_sw,fm_sw,'filled')
hline=refline(1,0);
hline.Color ='k';
hline.LineStyle='--';
hold on
plot(ff_sw,yCalc1)
ylabel('female to male');
xlabel('female to female');
ylim([0 1])
xlim([0 1])
hold off
[rf,pf]=corr([ff_sw; fm_sw]');
%functions

function f = smallangle(n)
    if ~isnan(n)
        
        if size(n,1)>1
            f = zeros(size(n,1),1);
            for val = 1:size(n,1)
                    if n(val) > pi
                        f(val) = (pi*2) -n(val);
                    else
                        f(val) = n(val);
                    end
            end
        else
            if n > deg2rad(180)
                f = (pi*2) - n;
            else
                f = n;
            end
        end
    else
        f = NaN();
        return
    end
end

function [f,gap_ext] = blocker(n,gap_closer) %put in matrix to translate to blocks; second argument tells the function to close gaps or not value is frames

    f = {}; %make cell array

    gap_ext = n;

    for fly = 1:size(n,2)
         
       blocks = zeros(1,4);
       start = NaN();
       for row =   2: size(n,1)
  
           if n(row,fly) >0 %is there interaction?
                  if n(row-1,fly) ~= n(row,fly) %is interaction new?
                     
                      if ~isnan(start) %was interaction with other fly before? i.e. is start defined?
                            blocks = [blocks; start row-1 n(row-1,fly) (row-1)-start+1]; %save range of last interaction
                           
                      
                      end
                      
                     start = row; %save start of current interaction
                     
                     if row == size(n,1) %if last one is singleton
                         blocks = [blocks; start row n(row,fly) 1];
                     end
                  else
                      if n(row-1,fly) == n(row,fly) && row == 2
                         start = row-1; 
                      else %if last one is interaction, save it
                          if row == size(n,1) && n(row,fly)>0
                              
                              blocks = [blocks; start row n(row,fly) (row)-start+1]; 
                          end
                      end
                  end
                  
           else
               if ~isnan(start) %if n(row,fly) is 0 but start is defined block needs to be saved
                  blocks = [blocks; start row-1 n(row-1,fly) (row-1)-start+1]; %save range of last interaction
                  start = NaN();
               end
           end
           
           
       end
        
       
       if gap_closer
          bl = 2;
          array_dim = size(blocks,1);
          while bl< array_dim
              if blocks(bl+1,1) - blocks(bl,2)<gap_closer && blocks(bl+1,3) == blocks(bl,3)%check if end of last block and begin of next block is less than gap_closer frames away from each other.
                  gap_ext(blocks(bl,2):blocks(bl+1,2),fly)=blocks(bl,3); %closes gap in map array
                  blocks(bl,2) = blocks(bl+1,2); %extends block length
                  blocks(bl+1,:) = []; %deletes obsolete row that was merged into previous block
                  blocks(bl,4) = blocks(bl,2)-blocks(bl,1)+1; %recalculates length of merged block
                  
                  
                  array_dim = array_dim-1;
              else
                  bl=bl+1;
              end
             
              
          end
       end
       
       f{fly,1} = blocks(2:end,:);
    end

end

function block_printer(n,numflies,adj,col,max_time,sex,min_lens)

    for l = 1:size(n,1)
        p = n{l};
        figure(1+l);
        yt = 1:numflies;
        yticks([yt]);
        ylabs = strcat('fly',num2str(yt'));
        yticklabels({ylabs});
        ylim([0 numflies]);
        for b = 1:size(p,1)
            if p(b,4)>=min_lens
                y=p(b,3)-adj;
                line([p(b,1),p(b,2)],[y,y],'LineWidth',5,'Color',col);
                y=l;
                col2=0;
                if(sex(p(b,3)))==1
                    col2='#FF9900';
                else
                    col2='#00AA00';
                end
                line([p(b,1),p(b,2)],[y,y],'LineWidth',5,'Color',col2);
            end
            
        end
        xlim([0 max_time]);
    end
end

function f = angle_difference(n,m)

    if sign(n) == sign(m)                

        f = abs(n - m) ;
    else
        f = smallangle(abs(n) + abs(m));
    end


end

function signals = detect_peaks(n,lag,threshold,influence)
   
    %n = smoothdata(n,'sgolay','SmoothingFactor', 0.03);
    mean_win_size=2000; %20000 for long videos - for 10 min videos not possible
    mean_win=mean_win_size;
    mean_all = nanmean(n(1:mean_win));
    signals = zeros(size(n,1),1);
    filteredY = n(1:lag);
    avgFilter(lag,1) = nanmean(n(1:lag));
    stdFilter(lag,1) = nanstd(n(1:lag));
    artifact = single.empty();

    for i=lag+1 : size(n,1)
        if(i>mean_win)
            if mean_win+mean_win_size>size(n,1)
               mean_win_size = size(n,1)-mean_win; 
            end
           mean_all = nanmean(n(mean_win+1:mean_win+mean_win_size));
           mean_win=mean_win+mean_win_size;
        end
       if abs(n(i) - avgFilter(i-1)) > threshold * stdFilter(i-1)
          if n(i) > avgFilter(i-1) && n(i) >1.5*avgFilter(i-1) && n(i)>mean_all
              signals(i) = 1;
              
              min_sig_win=50; %make sure not to run out of the array bounds here
              if(i-min_sig_win<1)
                 min_sig_win =i-1; 
              end
              
              if filteredY(i-1) > mean_all*2.5 && sum(signals(i-min_sig_win:i))>15
                  filteredY(i) = filteredY(i-1)*0.5;
              else
                filteredY(i) = influence * n(i) +(1-influence) * filteredY(i-1);
              end
              
              if n(i)>=1 %remove artefact peaks and adjust the signal matrix
                  if i< size(n,1)
                    artifact = [artifact i];
                  end
              end
          else
              filteredY(i) = n(i);
          end
          
       else
           signals(i) =0;
           filteredY(i) = n(i);
           
       end

       avgFilter(i) = nanmean(filteredY(i-lag:i));
       stdFilter(i) = nanstd(filteredY(i-lag:i) );
    end
    
    %remove artifact peaks
    if ~isempty(artifact)
        for err = 1:length(artifact)
           
            j=artifact(err);
            k=j;


            beginSignal=0;
            endSignal=0;
           
            
            while n(j)>n(j-1) %go backwards 'downhill' peak

                j=j-1;
                if(j==1)
                    
                   break 
                end
            end
            beginSignal = signals(j); %what is the state of Signal before artifact peak begins?

            while n(k)>n(k+1) %go forwards 'downhill' peak

                if(k+1<size(n,1))
                    k=k+1;
                else
                    break
                end

            end
            endSignal=signals(k); %what is the state of Signal after artifact peak ends

            %if endSignal and beginSignal are both ==1\, nothing needs to be
            %changed
            if endSignal ~= beginSignal %if either end or start signal are 0, signal at artifact peak needs to be changed to 0
                signals(j:k)=0;
            elseif endSignal == 0 && beginSignal ==0
                signals(j:k)=0;

            end
            
        end
    end
end

function map = conflict_resolver(n,win_size,dist) %win size to allow flexibility among behavioral motifs. solving mating ambiguities would benefit from a larger window than approaching or interacting
        %n is a cell array of numflies x 1 dimensions that contains a
        %frames x numflies array for each cell (i.e. numflies number of
        %cells)
        %what needs to happen is that whenever there is conflict in
        %approach, interaction, and mating calls these need to be solved
        %correctly. I'll start with the mating calls, because those are
        %easiest
        p=0; %print switch
        
        num_ind = size(n,1);
        
        map=zeros(size(n{1,1},1),num_ind);
        
        for i = 1:num_ind      
                              
                   curr = n{i,1}; %get target array from cell
                   
                   for step = 1:win_size:size(curr,1)
                                        
                      
                       to=step+win_size-1;
                       
                       if to>size(n{1,1},1)
                          to= size(n{1,1},1);
                       end
     if p                  
     fprintf('now at step %d to %d\n',step,to);
     end
                       
                       presence_tester=zeros(num_ind,1);
                       
                       for within = 1:num_ind
                           if sum(curr(step:to,within))
                                presence_tester(within,1) = 1; %indicates that a target was scored in the respective window
                           end
                       
                       end
                       
                       if sum(presence_tester) > 1 %more than 1 target individuals were scored for behavior. might need correction.
         
                           if p
                           fprintf('\tpotential conflict detected');
                           end
                           %now need to check if overlap
                           % 1111100000
                           % 0020022222
                           %or if successive
                           % 1111100000
                           % 0000022222
                           
                           %get a few values to decide what to keep and
                           %what to correct
                           
                           
                           targets = find(presence_tester);
                           if p
        fprintf(' - targets:')
           fprintf('%d ',targets);
             fprintf('\n');
                           end
                           num_calls=zeros(num_ind,1); % to identify number calls precedence
                           distances=num_calls; %distance precedence
                           past=num_calls; %past precedence

                           %translate matrix to blocks to test if there is
                           %overlap - gap closer is set to win size to get
                           %max span. This will lead to single block per
                           %individual.
                           [a_blocks,~] = blocker(curr(step:to,targets),win_size); 
                           
                           overlap=[];
                           
                           %join blocks saved in cell array in a
                           %concatenated array
                           total_blocks=0;
                           for cll =1:size(a_blocks,1)
                            total_blocks=total_blocks + size(a_blocks{cll,1},1);
                           end
                           action_blocks=zeros(total_blocks,4);
                           
                           las=1;
                           for ab = 1: size(a_blocks,1)
                                how_many=size(a_blocks{ab,1},1);
                                action_blocks(las:las+how_many-1,:) = a_blocks{ab,1};
                                las = how_many+las;
                           end
                         
                           for bl = 1:size(action_blocks,1)
                               
                               for bl2=bl+1:size(action_blocks,1)
                                  
                                   if action_blocks(bl,1) < action_blocks(bl2,1) && action_blocks(bl,2) < action_blocks(bl2,1)
                                        %no overlap
                                   elseif action_blocks(bl,1) > action_blocks(bl2,2)
                                        %no overlap
                                   else
                                       %overlap
                                       
                                       overlap = [overlap targets(bl) targets(bl2)];
                                   end
                               end
                             
                           end
                           
                           targets = unique(overlap); %only work on individuals that overlap
                           
                           if size(targets,1)
                               if p
       fprintf('\toverlap detected\n')
                               end
                               
                               for who = targets
                                   
                                   if p
            fprintf('\tworking on %d\n\t this is the sum %d\n\n',who,sum(curr(step:to,who)))
                                   end
                                   
                                  num_calls(who,1) = sum(curr(step:to,who))/who; %number of calls
                                  distances(who,1) = nanmean(dist{i,1}(step:to,who)); %mean distance in window for each individual

                                  from=step-win_size;

                                  if from<1
                                      from=1;
                                  end

                                  past(who,1) = sum(curr(from:step-1))/who; %checks who's the main interactor in past
                               end

                               [~,num_precedence]=sort(num_calls,'descend');
                               [~,dist_precedence]=sort(distances,'descend');
                               [~,past_precedence]=sort(past,'descend');
                               % figure out precedence hierarchy.
                               winner = 0;
                               if past_precedence(1) == num_precedence(1) && past_precedence(1) == dist_precedence(1)
                                   winner = past_precedence(1);

                               else
                                   % past_precedence wins over all others but only if num_precedence of some other individual is not much higher

                                   perc_calls=num_calls/max(num_calls);


                                   if perc_calls(past_precedence(1))< 0.75 % if individual with most past precedence is less then 75% of individual with maximum calls

                                       if num_calls(num_precedence(1))/sum(num_calls) < 2/3
                                           %if the fly with highest precedence
                                           %in number of frames interacting
                                           %makes up less than 2/3 of all
                                           %interactions of all individuals
                                           %combined

                                           if num_precedence(1) == dist_precedence(1) && num_calls(num_precedence(1))/sum(num_calls)>=0.5

                                               winner = dist_precedence(1);

                                           else

                                               sub_winner=0;
                                               for pr = 1:length(num_precedence)

                                                   if ~sub_winner
                                                      if num_calls(dist_precedence(pr))/sum(num_calls)>1/3

                                                          sub_winner = dist_precedence(pr);
                                                      end
                                                   end


                                               end

                                               if sub_winner
                                                    winner = sub_winner;

                                               else
                                                   winner = num_precedence(1);
                                               end

                                           end

                                       else

                                           winner = num_precedence(1);


                                       end

                                   else
                                       winner = past_precedence(1);

                                   end
                               end
                                if p
            fprintf('\tinteractor: %d\twinner: %d\n',i,winner)
                                end

                               map(step:to,i) = curr(step:to,winner);

                               %mp2=mating_cell{6,1}(25201:37800,2)
                               %%% this will also discard individuals that
                               %%% interacted before an overlap of individuals
                               %%% other than the first is detected. selection
                               %%% of good window size is key
                           else
                               if p
                               fprintf ('\tno overlap detected - all good\n')
                               end
                               t = find(presence_tester);
                               for put = 1:size(action_blocks,1)
                                   for tar = 1:t
                                       if action_blocks(put,3)==tar
                                        map(action_blocks(put,1):action_blocks(put,2))=curr(action_blocks(put,1):action_blocks(put,2),tar);
                                       end
                                   end
                               end
                           end
                          
                       else
                           
                           if sum(presence_tester)==1                                                     
                               atarget=find(presence_tester);
                               if p
             fprintf('\ntarget found: %d\n',atarget);
                               end
                                map(step:to,i)=curr(step:to,atarget);
                           
                           else
                              %nothing to do here, because no interactions, i.e. all 0 
                           end
                       end
                       
                   end
                    
                    
              
        end
end