function nii_mm2obj (Xmm,Ymm,Zmm, Radius, fnm, repo)
%Create mesh-based sphere at desired coordinates
% Xmm: Positon (Left-Right)
% Ymm: Position (Posterior-Anterior)
% Zmm: Position (Inferior-Superior)
% Radius: size of sphere (in mm)
% fnm: filename for sphere, e.g. 'test.obj'
% repo: (optional) reposition to surface of mesh. This refers to obj-format 
%        mesh: coordinates will be shifted to nearest vertex in this mesh
% mesh
% surface:
%Notes
% Chris Rorden, BSD 2-Clause license, 7/2016
%Examples
% nii_mm2obj(-51, -75, 0, 8, 'sphr')
% nii_mm2obj(-51, -75, 0, 8, 'sphr', 'T1.obj')
%

if exist('repo','var')
   [~, v] = readObjSub(repo);
    coord = [Xmm Ymm Zmm];
    dx = sqrt(sum(bsxfun(@minus,v,coord).^2,2)); %distance of each vertex
    [~,i] = min(dx); %index(es) of minimum
    Xmm = v(i(1),1);
    Ymm = v(i(1),2);
    Zmm = v(i(1),3);
    fprintf('Nearest vertex to %gx%gx%g is %gx%gx%g (%g)\n', coord(1), coord(2), coord(3), Xmm,Ymm,Zmm, min(dx));
end

vox=9;
thresh = 0.0;
[X,Y,Z]=ndgrid(linspace(-1,1,vox),linspace(-1,1,vox),linspace(-1,1,vox));
F = sqrt(X.^2 + Y.^2 + Z.^2);
F = 1.0 - F; %unit sphere distance from surface: center is 1, surface is 0, negative is outside
FV = isosurface(X, Y, Z, F, thresh); 
%FV = reducepatch(FV,0.2);%not needed:  vox=7 for coarse, vox=49 for fine 
if isempty(FV.vertices)
   fprintf('%s failed\n', mfilename)
    return;
end
%next lines add vertex color: shading based on position of vertex
clr = FV.vertices(:,1);
range = max(clr) - min(clr);
if range ~= 0 %normalize for range 0 (black) to 1 (white)
    FV.vertexColors = (clr - min(clr)) / range; %save colors as Scalar not RGB
end 
FV.faces = flipdim(FV.faces,2); %correct triangle winding
FV.vertices = FV.vertices * Radius;
FV.vertices(:,1) = FV.vertices(:,1) + Xmm;
FV.vertices(:,2) = FV.vertices(:,2) + Ymm;
FV.vertices(:,3) = FV.vertices(:,3) + Zmm;

[p,n] = fileparts(fnm);
writeObjSub(FV.vertices,FV.faces, fullfile(p,[n,'.obj']) );
%end nii_mm2obj()

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

function [faces, vertices] = readObjSub(fileName)
%function [faces, vertices] = readObj(fileName)
%inputs:
%	fileName: the Wavefront Object file to open
%outputs:
%	faces: face matrix where cols are xyz and each row is face
%	vertices: vertices matrix where cols are xyz and each row a vertix
%notes
% https://en.wikipedia.org/wiki/Wavefront_.obj_file

fid=fopen(fileName);
%first pass : count items
num_v=0;
num_f=0;
tline = fgets(fid);
while ischar(tline)
    A=strread(tline,'%s','delimiter',' ');
    if length(A) > 1
        if strcmpi(A(1),'v'), num_v = num_v + 1; end;
        if strcmpi(A(1),'f'), num_f = num_f + 1; end;
        %
    end
    tline = fgets(fid);
end
fclose(fid);
if (num_f < 1) || (num_v < 3)
    faces = [];
    vertices =[];
    fprintf('Unable to read this file');
    return; 
end;
fprintf('Obj file has %d vertices and %d faces\n', num_v, num_f);
%2nd pass : read items
fid=fopen(fileName);
num_v=0;
num_f=0;
vertices= zeros([num_v, 3]);
faces= zeros([num_f, 3]);
tline = fgets(fid);
while ischar(tline)
    A=strread(tline,'%s','delimiter',' ');
    if length(A) > 1
        if strcmpi(A(1),'v'), 
            num_v = num_v + 1;
            vertices(num_v, 1) = str2double(char(A(2)));
            vertices(num_v, 2) = str2double(char(A(3)));
            vertices(num_v, 3) = str2double(char(A(4)));
        end;
        if strcmpi(A(1),'f')
            num_f = num_f + 1;
            % we need to handle f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
            A2=strread(char(A(2)),'%s','delimiter','/');
            A3=strread(char(A(3)),'%s','delimiter','/');
            A4=strread(char(A(4)),'%s','delimiter','/');
            faces(num_f, 1) = str2double(char(A2(1)));
            faces(num_f, 2) = str2double(char(A3(1)));
            faces(num_f, 3) = str2double(char(A4(1)));
        end;
    end;
    tline = fgets(fid);
end
fclose(fid);
%end readObj()
