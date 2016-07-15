function nii_makeStatMask()
%188|PSIG_L|posterior inferior temporal gyrus left|1
%186|PSMG_L|posterior middle temporal gyrus left|1
%184|PSTG_L|posterior superior temporal gyrus left|1
%43|ITG_L|inferior temporal gyrus left|1
%39|MTG_L|middle temporal gyrus left|1
%35|STG_L|superior temporal gyrus left|1
%31|AG_L|angular gyrus left|1
%29|SMG_L|supramarginal gyrus left|1

roi = [29, 31, 35, 39, 43, 184, 186, 188];
fnm ='jhu.nii';
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
imgx = zeros(size(img));
str ='';
for i = 1 : numel(roi)
	imgx(img == roi(i)) = 1;
	str = [str, ' ', num2str(roi(i))];
end;
hdr.fname = 'mask.nii';
spm_write_vol(hdr,imgx);
nii_reslice_target('mask.nii',[],'statmask.nii',0);
fprintf('"rmask.nii" composed of regions %s\n', str);
fprintf('Rename "rmask.nii" as "statmask.nii" to alter region of interest\n');
