% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% computes Tracking Index of visual stimulus bound between -1 and 1
%
% function used by 'Brand_et_al_26_functional_imaging_postProcessing.m'
% rename this file to 'computeTI.m' for it to work
% function modified after Hindmarsh Sten et al. 2021 (Hindmarsh Sten, T.,
% Li, R., Otopalik, A., & Ruta, V. (2021). Sexual arousal gates visual processing during Drosophila courtship. Nature, 595(7868), 549-553.)
% cite Hindmarsh Sten et al. 2021 if you use this script or a modified version


function [R,vigor,fidelity] = computeTI(hDiff,femaleAngle,bTime)

%parameters
trialLength = 180; %bin size 
hemiLength = trialLength/2; %half-length
trialStarts = 1:length(bTime)-1; %frameIDs
ctr = 1;

%compute tracking index
for i = trialStarts
    
    %if not at edges of trial
    if i > hemiLength && i < (length(bTime)-hemiLength)
        x = hDiff((i-hemiLength):(i+hemiLength)); %male turning
        f = femaleAngle((i-hemiLength):(i+hemiLength));%stim pos
        t = bTime((i-hemiLength):(i+hemiLength)); %time stamps
        cTime(ctr) = bTime(i); %time 
        mat = corrcoef(x,f); %correlation matrix between stim and male turning
        fid = mat(2); %store correlation coefficient
        
    % if at the very start of trial (too few n), set to zero 
    elseif i <= hemiLength
        x = hDiff(1:i);
        f = femaleAngle(1:i);
        t = bTime(1:i);
        cTime(ctr) = bTime(i);
        fid = 0; 
    
    % if at the very end of trial (too few n), set to zero 
    else
        x = hDiff(i:end);
        f = femaleAngle(i:end-1);
        t = bTime(i:end-1);
        cTime(ctr) = bTime(i);
        fid = 0; 
    end
    
    fidelity(ctr) = fid; %correlation
    vigor(ctr) = sum(x.*sign(f)); %net turning in direction of target
    ctr = ctr+1;
end

R = vigor./max(vigor); %normalize vigor across all periods
R = R.*(fidelity); %tracking index