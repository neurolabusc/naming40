function nii_block_2sess_jhu (t1, fmri)
%sample analysis of a block design with two sessions of data

if ~exist('t1','var')  %no files specified
 t1 = spm_select(1,'image','Select T1 scan');
end
if isempty(t1) || (size(t1,1) ~= 1)  %no files specified
 error('No T1 image');
end
if ~exist('fmri','var')  %no files specified
 fmri = spm_select(inf,'image','Select fMRI scan(s)');
end
if isempty(fmri) %no files specified
 error('No fmri image');
end
if (size(fmri,1) < 1) ||  (size(fmri,1) > 2)
    error('Select either one or two fMRI sessions');
end

p.fmriname = fmri; %strvcat('fmriblocks009a.nii', 'fmriblocks009b.nii'); %#ok<REMFF1>
p.t1name = t1; %'t1.nii'
p.TRsec = 10.0; %repeat time off 10 seconds
p.slice_order = -1; %skip slice-time corrections
p.phase =  ''; %phase image from fieldmap
p.magn = ''; %magnitude image from fieldmap
%statistical information (optional: if not provided not statitics)
p.names{1} = 'STIM';
p.names{2} = 'ABS';
p.forceTemporalDeriv = true;
%onset timing
p.onsets{1,1} = [16.734 26.031 40.609 51.641 70.141 88.625 109.469 128.516 138.75 169.813 186.813 201.141 208.859 219.875 250.672 259.875 266.359 278.953 307.234 328.422 337.219 349.531 358.922 367.844 377.547 388.234 401.109 410.891 421 440.406 451.469 468.563 477.625 506.578 520.594 541.813 549.75 581.125 587.406 597.563 606.344 627.672 650.031 677.734 700.969 718.578 726.219 739.5 747.672 758.453 791.125 796.953 808.609 821.313 831.313 839.078 848.781 860.188 886.781 906.734 926.219 947.406 956.938 970.391 976.734 989.875 1000.859 1021.453 1030.141 1049.969 1057.281 1068.656 1086.891 1100.781 1106.641 1117.781 1128.344 1135.938 1158.031 1170.75];
p.onsets{1,2} = [11.547 57.813 80.797 96.266 119.016 147.891 158 180.672 228.016 237.36 287.157 301.813 320.891 426.75 456.75 491.704 496.672 530.641 556.344 569.829 616.047 641.954 661.875 669.235 689.547 711.141 767.954 780.532 870.141 877.016 901.532 917.047 941.844 1008.797 1041.329 1078.016 1148.079 1177.954 1187.422 1200.704];

%timings are relative to middle of volume (FSL timing for Sparse, adjust to actual times)
p.onsets{1,1} = p.onsets{1,1} - (p.TRsec/2);
p.onsets{1,2} = p.onsets{1,2} - (p.TRsec/2);
if size(fmri,1) == 2
    %onsets for 2nd session, onsets{session, cond}
    p.onsets{2,1} = p.onsets{1,1};
    p.onsets{2,2} = p.onsets{1,2};
end


%kDur = 3; %each trail is 167ms target plus 1600ms fixation, 16 trials per block
p.duration{1} = 3; %duration 1 for events, longer for blocks - we have used 3sec for sparse blocks
p.mocoRegress = true; %should motion parameters be included in statistics?
%save settings
pth = fileparts(p.t1name);
cd(pth);
save('naming40.mat', '-struct', 'p');
%run the analysis
nii_batch12_jhu(p);

