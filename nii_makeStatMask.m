function nii_makeStatMask()
%188|PSIG_L|posterior inferior temporal gyrus left|1
%186|PSMG_L|posterior middle temporal gyrus left|1
%184|PSTG_L|posterior superior temporal gyrus left|1
%39|MTG_L|middle temporal gyrus left|1
%35|STG_L|superior temporal gyrus left|1
%31|AG_L|angular gyrus left|1
%29|SMG_L|supramarginal gyrus left|1
%21July2016No longer included:
% 43|ITG_L|inferior temporal gyrus left|1



roi = [29, 31, 35, 39, 184, 186, 188];
%roi = [29, 31, 35, 39, 43, 184, 186, 188];

pth = fileparts(which(mfilename));
fnm = fullfile(pth, 'jhu.nii');
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
imgx = zeros(size(img));
str ='';
for i = 1 : numel(roi)
	imgx(img == roi(i)) = 1;
	str = [str, ' ', num2str(roi(i))];
end;
tempfnm = fullfile(pth, 'temp.nii');
maskfnm = fullfile(pth, 'statmask.nii');
oldfnm = fullfile(pth, 'statmask_old.nii');
hdr.fname = tempfnm;
spm_write_vol(hdr,imgx);
nii_reslice_target(tempfnm,[],maskfnm,0);
delete(tempfnm);
tempfnm = fullfile(pth, 'rtemp.nii'); %resliced
copyfile(maskfnm,oldfnm);
movefile(tempfnm,maskfnm);
fprintf('updated %s\n', maskfnm);
