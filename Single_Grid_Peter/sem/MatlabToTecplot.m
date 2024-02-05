function []=MatlabToTecplot(style,filename,title,variablenames,meshsize,data,dimension,corners)
%
% MatlabToTecplot(style,filename,title,variablenames,meshsize,data,dimension,corners)
%
% This function creates ASCII filename.dat files as data files for
% TecPlot360.
%
% style : This function can produce tecplot data-files in two different
% styles, ordered and finite element style, indicated by either
%             style = 'IJK'
%             style = 'FE'
%
% filename: is the name of the .dat file
% title   : will be read Tecplot
% variablenames: names of the variables, given as a string, where the
% variables are double quoted and separated by spaces, e.g.
%             variablenames = '"A" "B" "C"'
% meshsize: is an array that either indicated the number of points in I,J
% and K direction, or indicates the number of nodes and number of cells
% data: Set of data columns
% dimension: dimension of the problem, i.e. 1D, 2D, 3D
% corners: (only needed for style='FE') indicates the order of nodes of a
% cell


% creation of filename
filename = strcat(filename,'.dat');

if exist(filename,'file')
    disp(['existing file ' filename ' is replaced by a new file'])
    if ispc
        system(['del ' filename]);
    elseif isunix
        system(['rm ' filename]);
    end
end

% Open file
fid = fopen(filename,'w');

% writing title and variablenames
fprintf(fid,['TITLE = "' title '"\n']);
fprintf(fid,['VARIABLES = ' variablenames '\n']);

% Writing style info
switch style
    
    case 'IJK'
        
if dimension==1
fprintf(fid,['ZONE ZONETYPE=ORDERED, I=' num2str(meshsize(1)) ', DATAPACKING=POINT\n']);
elseif dimension==2
fprintf(fid,['ZONE ZONETYPE=ORDERED, I=' num2str(meshsize(1)) ', J=' num2str(meshsize(2)) ', DATAPACKING=POINT\n']);
elseif dimension==3
fprintf(fid,['ZONE ZONETYPE=ORDERED, I=' num2str(meshsize(1)) ', J=' num2str(meshsize(2)) ', K=' num2str(meshsize(3)) ', DATAPACKING=POINT\n']);
end

    case 'FE'

if dimension==1
fprintf(fid,['ZONE N=' num2str(meshsize(1)) ' E=' num2str(meshsize(2)) ', DATAPACKING=POINT, ZONETYPE=FELINESEG\n']);
elseif dimension==2
fprintf(fid,['ZONE N=' num2str(meshsize(1)) ' E=' num2str(meshsize(2)) ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n']);
elseif dimension==3
fprintf(fid,['ZONE N=' num2str(meshsize(1)) ' E=' num2str(meshsize(2)) ', DATAPACKING=POINT, ZONETYPE=FEBRICK\n']);
end

end

% writing data
str = '%15.8f';
for i=1:size(data,2)-1
    str = strcat(str,'\t%15.8f');
end
str = strcat(str,'\n');

fprintf(fid,str,data');
fprintf(fid,'\n');

% writing corners
if strcmp(style,'FE') && exist('corners','var')
str = '%8.u';
for i=1:size(corners,2)-1
    str = strcat(str,'%8.u');
end
str = strcat(str,'\n');
fprintf(fid,str,corners');
end

%closing file
fclose(fid);