function naming40(imgDir)
%process sparse fMRI data from JHU
% imgDir : (optional) folder containing *MPrage*PAR and *fMRIsparse*PAR images
%Examples
% naming40
% naming40(pwd)

fprintf('%s version 21July2016\n', mfilename);
checkForUpdateSub(fileparts(mfilename('fullpath')));
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running
if ~exist('imgDir','var')
    title = 'Select folder with images (fMRI*.nii and T1*.nii or raw PAR/REC files)';
    fprintf('%s\n', title);
    imgDir = uigetdir(pwd, title);
    %imgDir = pwd; %666
end;
cd(imgDir);
if ~exist(imgDir,'file'), fprintf('Unable to find %s\n', imgDir); return; end;
matname = fullfile(imgDir, 'naming40.mat');
if ~exist(matname, 'file') %must compute analyses
   t1 = findImgsSub(imgDir, 'T1');
   fmri = findImgsSub(imgDir, 'fmri');
   if isempty(t1) || isempty(fmri)
        par2niiSub;
        t1 = findImgsSub(imgDir, 'T1');
        fmri = findImgsSub(imgDir, 'fmri');
        if isempty(t1) || isempty(fmri)
           return;
        end
   end;
   nii_block_2sess_jhu (t1, fmri)
else
    fprintf(' Skipping analyses: to re-compute delete %s\n', matname);
end
m = load(matname);
%if false %666
if  isfield(m, 'scriptname') && exist(m.scriptname, 'file');
    fprintf(' Skipping mesh generation: to re-compute delete %s\n', m.scriptname);
else
    makeMeshesSub(matname, false); %in warped space (optional)
    makeMeshesSub(matname, true); %in native space

end
if ~isfield(m, 'surficename') || ~exist(m.surficename, 'file');
    findSurficeSub(matname);
    m = load(matname);
end
sys = sprintf('%s -s "%s"', m.surficename, m.scriptname);
if isunix
    system([sys, ' &']);
else
    disp('For windows, perhaps asynchronous: http://stackoverflow.com/questions/18111211/use-matlab-to-call-system-function-asynchronously-in-background');
    system(sys);
end
%end nii_naming40()

function par2niiSub
t1Key = '*MPrage*PAR';
fmriKey = '*fMRIsparse*PAR';
t1Dir = dir(t1Key);
fmriDir = dir(fmriKey);
if isempty(t1Dir), fprintf('Unable to find file matching %s\n', t1Key); return; end;
if isempty(fmriDir), fprintf('Unable to find file matching %s\n', fmriKey); return; end;
t1 = fullfile(pwd, t1Dir(1).name);
fmri = fullfile(pwd, fmriDir(1).name);
if numel(t1Dir) > 1, warning('Multiple files: will use %s', t1); end;
if numel(fmriDir) > 1, warning('Multiple files: will use %s', fmri); end;
dcm2niix = findDcm2NiiSub;
sys = sprintf('%s -f fmri "%s"', dcm2niix, fmri);
system(sys);
sys = sprintf('%s -f t1 "%s"', dcm2niix, t1);
system(sys);
%end par2niiSub();

function dcm2niix = findDcm2NiiSub
p = fileparts(which(mfilename));
if ismac
    dcm2niix = fullfile(p, 'dcm2niix');
elseif isunix
    dcm2niix = fullfile(p, 'dcm2niixLX');
else
    dcm2niix = fullfile(p, 'dcm2niix.exe');
end
if ~exist(dcm2niix, 'file'), error('Unable to find %s', dcm2niix); end;
%findSurficeSub()


function findSurficeSub(matname)
p = fileparts(which(mfilename));
m = load(matname);
if ismac
    m.surficename = fullfile(p, '/surfice.app/Contents/MacOS/surfice');
elseif isunix
    m.surficename = fullfile(p, 'surfice');
else
    m.surficename = fullfile(p, 'surfice.exe');
end
if ~exist(m.surficename, 'file'), error('Unable to find %s', m.surficename); end;
save(matname, '-struct', 'm');
%findSurficeSub()

function makeMeshesSub(matname, nativeSpace)
%Make surface of scalp and brain, show sphere at peak activation on brain and nearest scalp
m = load(matname);
statmask = fullfile(fileparts(which(mfilename)),'statmask.nii');
if ~exist(statmask, 'file'), error('Unable to find %s', statmask); end;
%find peak statistics
[p,n] = fileparts(m.fmriname);
statname = fullfile(p,n,'spmT_0001.nii');
if ~exist(statname,'file'), fprintf('Unable to find %s\n', statname); return; end;
[~, mask] = readImgSub(statmask);
[hdr, img] =  readImgSub(statname);
img(mask== 0) = 0;
%next lines unique to Native-space version
prefix = 'w'; %use normalized ('w'arped space)
if (exist('nativeSpace','var')) && (nativeSpace)
    prefix = ''; %use native images (not 'w'arped space)
    %save masked image in warped space
    maskstatname = fullfile(p,n,'mskspmT_0001.nii');
    hdr.fname = maskstatname;
    spm_write_vol(hdr,img);
    [p,n] = fileparts(m.t1name);
    inv = fullfile(p,['iy_', n, '.nii']); %inverse
    mb{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {inv};
    mb{1}.spm.util.defs.comp{1}.inv.space = {inv};
    mb{1}.spm.util.defs.out{1}.push.fnames = {maskstatname};
    mb{1}.spm.util.defs.out{1}.push.weight = {''};
    mb{1}.spm.util.defs.out{1}.push.savedir.savepwd = 1;
    mb{1}.spm.util.defs.out{1}.push.fov.file = {inv};
    mb{1}.spm.util.defs.out{1}.push.preserve = 0;
    mb{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    %deform from warped mni to native space
    spm_jobman('run',mb);
    %the prefix 'w' is misleading, this image have been inverted to native space
    maskstatnamex = fullfile(p,'wmskspmT_0001.nii');
    maskstatname = fullfile(p,'invspmT_0001.nii');
    movefile(maskstatnamex, maskstatname);
    %load native masked image
    [hdr, img] =  readImgSub(maskstatname);
end
%return to common code:
[mx, i] = max(img(:));
[vx(1) vx(2) vx(3)] = ind2sub(size(img),i(1));
v2m = hdr.mat; %voxel2mm transform
mm= vx*v2m(1:3,1:3)' + v2m(1:3,4)';
fprintf('Maximum statistical intensity %g (%gx%gx%g)\n', mx, mm(1), mm(2), mm(3));
if mx <= 0, error('No positive values in mask!'); end;
%create meshes
[p,n] = fileparts(m.t1name);
bet = fullfile(p,[prefix, 'b', n, '.nii']);
if exist(bet, 'file')
    m.brainname = fullfile(p,[prefix 'brain.obj']);
    nii_nii2obj (bet, m.brainname);
end;
sumname = makeSumSub(m.t1name, nativeSpace);
m.scalpname = fullfile(p, [prefix 'scalp.obj']);
nii_nii2obj (sumname, m.scalpname);
m.peakname = fullfile(p,[prefix 'peak.obj']);
nii_mm2obj (mm(1), mm(2), mm(3), 8, m.peakname);
m.surfacename = fullfile(p,[prefix 'surface.obj']);
nii_mm2obj (mm(1), mm(2), mm(3), 8, m.surfacename, m.scalpname);
m.scriptname = fullfile(p,[prefix 'scalp.gls']);
save(matname, '-struct', 'm');
glscript = sprintf('BEGIN\n RESETDEFAULTS;\n  COLORBARVISIBLE(false);\n MESHLOAD(''%s'');\n OVERLAYLOAD(''%s'');\n OVERLAYLOAD(''%s'');\n OVERLAYLOAD(''%s'');\n SHADERXRAY(1.0, 0.2);\nEND.', ...
     m.scalpname, m.peakname, m.surfacename, m.brainname);
fileID = fopen(m.scriptname,'w');
fprintf(fileID,glscript);
fclose(fileID);
%end makeMeshesNativeSub()

function [hdr, img] = readImgSub(fnm)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
img(isnan(img)) = 0;
% end readImgSub()

function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function nameFiles=subImgSub(pathFolder)
nameFiles=subFileSub(pathFolder);
if isempty(nameFiles), return; end;
n = nameFiles; nameFiles = [];
for i = 1: numel(n)
    [~,~,x] = fileparts(char(deblank(n(i))));
    if ~strncmpi('.gz',x, 3) && ~strncmpi('.nii',x, 4), continue; end;
    nameFiles = [nameFiles; n(i)]; %#ok<AGROW>
end
%end subFileSub()

function fnm = findImgsSub( xDir, xKey)
fnm = [];
nameFiles = subImgSub(xDir);
if isempty(nameFiles), return; end;
n = 0;
for j = 1: numel(nameFiles)
        if strncmpi(xKey,nameFiles(j),numel(xKey))
           if n == 0, fnm = fullfile(xDir, char(nameFiles(j)) ); end;
           n = n + 1;
        end
end;
if n == 0, error('Unable to find %s*.nii in %s', xKey, xDir); end;
if n > 1, warning('Found more than one %s*.nii in %s: using %s\n', xKey, xDir, fnm); end;
%end findImgsSub()

function checkForUpdateSub(repoPath)
prevPath = pwd;
cd(repoPath);
if exist('.git','dir') %only check for updates if program was installed with "git clone"
    [s, r] = system('git fetch origin','-echo');
    if strfind(r,'fatal')
        warning('Unabe to check for updates. Network issue?');
        return;
    end
    [~, r] = system('git status','-echo');
    if strfind(r,'behind')
        if askToUpdate
            [~, r] = system('git pull','-echo');
            showRestartMsg
        end
    end
else %do nothing for now
    warning(sprintf('To enable updates run "!git clone git@github.com:neurolabusc/%s.git"',mfilename));
end
cd(prevPath);
%end checkForUpdateSub()

function sumname = makeSumSub(t1name, nativeSpace)
[pth nm ext] = spm_fileparts(t1name);
c1 = fullfile(pth,['c1',nm,ext]);
c2 = fullfile(pth,['c2',nm,ext]);
c3 = fullfile(pth,['c3',nm,ext]);
c4 = fullfile(pth,['c4',nm,ext]);
c5 = fullfile(pth,['c5',nm,ext]);
if ~exist(c1, 'file') || ~exist(c2, 'file') || ~exist(c3, 'file') ||  ~exist(c4, 'file') || ~exist(c5, 'file')
    error('Unable to find tissue maps %s..%s', c1, c5);
end
hdr = spm_vol(c1);
img = spm_read_vols(hdr);
hdr = spm_vol(c2);
img = img + spm_read_vols(hdr);
hdr = spm_vol(c3);
img = img + spm_read_vols(hdr);
hdr = spm_vol(c4);
img = img + spm_read_vols(hdr);
hdr = spm_vol(c5);
img = img + spm_read_vols(hdr);
img(isnan(img)) = 0;
%save sum image
sumname = fullfile(pth,'sum.nii');
hdr.fname = sumname;
spm_write_vol(hdr,img);
if exist('nativeSpace', 'var') && (nativeSpace), return; end;
newSegWriteSub(t1name, sumname);
sumname = fullfile(pth,'wsum.nii');
%end readImgSub()


function  newSegWriteSub(t1name, warpname, resliceMM, bb)
%reslice img using pre-existing new-segmentation deformation field
if isempty(warpname) || isempty(t1name), return; end;
if ~exist('resliceMM','var'), resliceMM = 1; end;
[pth,nam,ext, vol] = spm_fileparts(t1name);
defname = fullfile(pth,['y_' nam ext]);
if ~exist(defname,'file')
    defname = fullfile(pth,['y_e' nam ext]);
    if ~exist(defname,'file')
        defname2 = fullfile(pth,['y_' nam ext]);
        error('Unable to find new-segment deformation image %s or %s',defname2, defname);
    end
end
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {defname};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {warpname};
if ~exist('bb','var')
	matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-91 -126 -76; 91 96 108]; %loose
    %matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85]; %tight
else
	matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
end
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [resliceMM resliceMM resliceMM];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
spm_jobman('run',matlabbatch);
%end newsegwritesub()