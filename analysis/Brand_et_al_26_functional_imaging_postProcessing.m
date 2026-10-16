% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% computes basic tracking parameters from virmen output
% plots male tracking of visual stimulus
%
% uses computeTI function in Brand_et_al_26_computeTI.m
%
% script modified after Hindmarsh Sten et al. 2021 (Hindmarsh Sten, T., Li, R., Otopalik, A., & Ruta, V. (2021). Sexual arousal gates visual processing during Drosophila courtship. Nature, 595(7868), 549-553.)
% cite Hindmarsh Sten et al. 2021 if you use this script or a modified version


%requires CircStat2012a to work
addpath('/Applications/MATLAB_R2022a.app/toolbox/CircStat2012a/')
% organize data and get basic params

invariant=1; %1 if invatiant visual stimulus was used

bTime = expr(:,1); %time stamp
male.pos = expr(:,3:4);  %male position in ViRMEn coords
male.heading = expr(:,5); %male heading in ViRMEn coords
female.pos = expr(:,[6:7] ); %female position in ViRMEn coords [6:7] default

% Define the rotation angle in degrees and convert to radians

shift=0;
if shift==1
        x=expr(:,6);
        y=expr(:,7);
        
        angle_degrees = 0; %to the left: 75 degrees, to the right: 285
        angle_radians = deg2rad(angle_degrees);
        
        % Create the rotation matrix
        %R = [cos(angle_radians), -sin(angle_radians); sin(angle_radians), cos(angle_radians)];
        
        % Apply the rotation to each point
        %rotated_points = R * [x; y];
        
        female.pos = [x*cos(angle_radians)-y*sin(angle_radians), x*sin(angle_radians)+y*cos(angle_radians)];
end

male.distance = sqrt(sum((female.pos).^2,2)); %distance between male and female
male.vel = sqrt(sum(diff(male.pos).^2,2)).*(1./diff(bTime)).*3; %male velocity
female.vel = sqrt(sum(diff(female.pos).^2,2)); %female velocity
male.distTraveled=sqrt(sum(diff(male.pos).^2,2)).*3;
zCoors = female.pos - male.pos; %zero coordinate system w/ respect to spectator
female.absAng = atan2(zCoors(:,2),zCoors(:,1));   %angular location of female w/ resp. to male

male.normHeading = circ_dist(male.heading,wrapTo2Pi(female.absAng-pi/2)); %heading with respect to female (mirror inversion corrected)
female.viewPosition = [male.distance.*cos(male.normHeading) male.distance.*sin(male.normHeading)]; %pos of female in male's FOV

%compute male's trajectory in VR (had he not been stationary)
transform = @(deltaXY,heading) [cos(heading) -sin(heading); sin(heading) cos(heading)]*deltaXY';
deltaPos = diff(male.pos);
male.VRpos(1,:) = [0 0];
for i = 2:length(expr)
    trueShift = transform(deltaPos(i-1,:),male.heading(i));
    male.VRpos(i,:) = male.VRpos(i-1,:) + trueShift';
end 


%compute angle of object;
angleOL = atan2(expr(:,7),expr(:,6));
female.angle = atan2(female.pos(:,1),female.pos(:,2));

%figure out how long the phase of the stimulus is
female.angle = round(female.angle,5);
firstMotion = find(diff(female.angle)~=0,1);

if invariant
    offset=1;
else
    offset=0;
end
phaseLength = find(female.angle(firstMotion:end,1) == female.angle(firstMotion,1))+offset;%add +1 for invar stim
phaseLength = phaseLength(2)*3;

%find number of full phases completed
nPrePhases = floor(firstMotion./phaseLength);
nPostPhases = floor((length(expr)-firstMotion)./phaseLength);
totalPhases = nPrePhases+nPostPhases;
firstID = firstMotion-nPrePhases*phaseLength;
lastID = firstMotion + nPostPhases*phaseLength-1;

%compute ?Heading of male
hDiff = zeros(length(expr)-1,1);
male.angVel = hDiff;
headTMP=male.heading;
for i = 1:length(expr)-1
    hDiff(i) = circ_dist2(male.heading(i),male.heading(i+1));
    if(headTMP(i+1)>pi)
        headTMP(i+1)=pi-(headTMP(i+1)-pi);
    end
    if(headTMP(i)>pi)
        headTMP(i)=pi-(headTMP(i)-pi);
    end
    
end
   male.angVel=diff(headTMP);
   male.angVel=male.angVel*(length(male.angVel)/bTime(end)); %to get to rad per s for angular velocity
   for i=1:length(male.angVel)
       
        if(male.angVel(i)<0)
            male.angVel(i)=male.angVel(i)*-1;
        end
   end






%align parameters into sets of three stim cycles
meanPhase = reshape(hDiff(firstID:lastID),phaseLength,totalPhases)';
stimPhase = reshape(female.angle(firstID:lastID),phaseLength,totalPhases)';
timePhase = reshape(bTime(firstID:lastID),phaseLength,totalPhases)';
meanHeading = reshape(male.normHeading(firstID:lastID),phaseLength,totalPhases)';

%get tracking index
[r,vig,fidelity] = computeTI_unbiased_bound(hDiff,female.angle,bTime);

%plot male turning over the course of the trial
map=[1 0 0
    0.8 0 0
    0.6 0 0
    0.4 0 0
    0.2 0 0
    0 0 0
    0 0 0.2
    0 0 0.4
    0 0 0.6
    0 0 0.8
    0 0 1];
figure;
subplot(121);
imagesc(meanPhase); caxis([-0.1 0.1]);
ylabel('Iteration','FontSize',15,'FontName','Avenir Book');
xlabel('Time (frame, 59Hz)','FontSize',15,'FontName','Avenir Book');
set(gca,'FontSize',12,'FontName','Avenir Book')
box off; colormap(map) %parula

%plot tracking index
subplot(122);
plot(r,bTime(2:end),'b','LineWidth',1.5)
set(gca, 'YDir','reverse')
ylabel('Time (sec)','FontSize',15,'FontName','Avenir Book')
xlabel('Tracking Index','FontSize',15,'FontName','Avenir Book')
ylim([min(bTime) max(bTime)])
xlim([-1 1])
xline(0)
set(gca,'FontSize',12,'FontName','Avenir Book')
box off;

set(gcf,'Color','w','Position',[441 459 486 491])

ti=[r',fidelity',vig',bTime(2:end)];


res=[bTime(2:end), female.angle(2:end), hDiff, r', fidelity',vig',male.vel, male.angVel];
%% save tracking parameters
%writematrix(res,'file_name.csv')

%save('file_name.mat','ti')




