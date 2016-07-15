function outnm = nii_nii2obj (fnm, outnm, reduce)
%convert NIfTI image to mesh
% fnm : nifti image to meshify
% outname : (optional) name for output mesh 
% reduce : (optional) reduction ration, e.g. 0.2 will decimate 80% of triangles
%Examples
% nii_nii2obj('c1T1.nii');

if ~exist('reduce','var')  %no files specified
    reduce = 0.08; %reduce mesh complexity to 10%, smaller values for smaller files
end
if ~exist('fnm','var')  %no files specified
    fnm = spm_select(1,'image','Select scan to meshify');
end
if ~exist('outnm', 'var')
    [pth nm] = spm_fileparts(fnm);
    outnm = fullfile(pth, [nm '.obj']);
end;
%load image
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
img(isnan(img)) = 0;
%blur data
presmooth = img+0; %+0 forces new matrix
spm_smooth(presmooth,img,6,0); %smooth data
%find threshold
thresh = (0.5 * max(img(:))-min(img(:))) + min(img(:));
fprintf('%g threshold for %s\n', thresh, fnm);
img = floodSub(img, thresh, 2);
FV = isosurface(img,thresh);
FV = reducepatch(FV,reduce);
FV.vertices = FV.vertices(:,[2,1,3]); %isosurface swaps X/Y
vx = [ FV.vertices ones(size(FV.vertices,1),1)];
vx = mtimes(hdr.mat,vx')';
FV.vertices = vx(:,1:3);
FV.faces = fliplr(FV.faces); %reverse winding
%save results
writeObjSub(FV.vertices,FV.faces, outnm);
%end nii_nii2obj()

function writeObjSub(vertex,face,filename)
% --- save Face/Vertex data as WaveFront Object format file
%inputs:
%	vertex: vertices matrix where cols are xyz and each row a vertix
%	face: face matrix where cols are xyz and each row is face
%	fileName: the Wavefront Object file to create
%notes
% https://en.wikipedia.org/wiki/Wavefront_.obj_file
[nF nFd] =size(face);
[nV nVd] =size(vertex);
if (nF <1) || (nV <3 || (nFd ~=3) || (nVd ~=3)), warning('Problem with writeObj'); return; end; 
fid = fopen(filename, 'wt');
fprintf(fid, '# WaveFront Object format image created with MRIcroS\n');
fprintf(fid, 'v %.12g %.12g %.12g\n', vertex');
fprintf(fid, 'f %d %d %d\n', (face)');
fclose(fid);
%end writeObjSub()


function img = floodSub(img, thresh, smoothVox)
%similar to http://www.mathworks.com/matlabcentral/fileexchange/12184-floodfill3d

imgBin = (img < thresh);
imgBin = double(imgBin);      % In case a logical matrix comes in.
imgBin(1,:,:) = NaN; %pad LEFT 
imgBin(end,:) = NaN; %pad RIGHT
imgBin(:,1,:) = NaN; %pad ANTERIOR 
imgBin(:,end,:) = NaN; %pad POSTERIOR
imgBin(:,:,1) = NaN; %pad INFERIOR 
imgBin(:,:,end) = NaN; %pad SUPERIOR
imgBin = floodFill3DSub(imgBin, [2,2,2]);
%fill bubbles
mx = max(img(:));
img(isfinite(imgBin))= mx;
%next: optional - blur to feather edges - only useful if marching cubes uses subvoxel smoothing
maskIn = img + 0.0;
if exist('smoothVox', 'var')
    spm_smooth(maskIn, img, smoothVox); %feather the edges a lot: weak blur
end
%end floodSub()

%Copyrights for Matlab Central Code
% stlwriteSub Copyright (c) 2015, Sven Holcombe
% floodFill3DSub Copyright (c) 2006,  F Dinath
% 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [A] = floodFill3DSub(A, point)
% http://www.mathworks.com/matlabcentral/fileexchange/12184-floodfill3d
% [B] = FloodFill3D(A, slice);
% This program flood fills a 6-connected 3D region. The input matrix MUST
% be a binary image. The user will select a seed (point) in the matrix to
% initiate the flood fill. You must specify the matrix slice in which you
% wish to place the seed.
% 
% A = binary matrix
% slice = a chosen slice in the matrix where you wish to place the seed.
%
% Enjoy,
% F. Dinath
if A(point(1), point(2), point(3));
    A(point(1), point(2), point(3)) = NaN;
    a{1} = sub2ind(size(A), point(1), point(2), point(3));
    i = 1;
    while 1
        i = i+1;
        a{i} = []; %#ok<AGROW>
        [x, y, z] = ind2sub(size(A), a{i-1});
        ob = nonzeros((A(sub2ind(size(A), x, y, z-1)) == 1).*sub2ind(size(A), x, y, z-1));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y, z+1)) == 1).*sub2ind(size(A), x, y, z+1));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x-1, y, z)) == 1).*sub2ind(size(A), x-1, y, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x+1, y, z)) == 1).*sub2ind(size(A), x+1, y, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y-1, z)) == 1).*sub2ind(size(A), x, y-1, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        ob = nonzeros((A(sub2ind(size(A), x, y+1, z)) == 1).*sub2ind(size(A), x, y+1, z));
        A(ob) = NaN;
        a{i} = [a{i} ob']; %#ok<AGROW>
        if isempty(a{i});
            break;
        end
    end
end
%end floodFill3DSub()