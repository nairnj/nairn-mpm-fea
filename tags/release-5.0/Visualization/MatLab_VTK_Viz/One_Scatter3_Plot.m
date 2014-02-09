% MatLab Script: One_Scatter3_Plot.m
%
% Date: 15 Feb 2011
% Author: John Nairn
%
% This script reads all data in a vtk file. You can select on
% quantity to plotted in a scatter plot. If the selected quantity
% is a tensor or vector, you can select which component to plot
%
% Tensors: 1 to 9 are xx,xy,xz,yx,yy,yz,zx,zy,zz
% Tensors: 1 to 3 are x, y, and z

% User select a file
[fileName,fldrPath,FilterIndex] = uigetfile('*.vtk','Select VTK File');
if(fileName==0)
    return
end

% read all data
[VTKheader,VTKdata] = ReadVTKFile(fldrPath,fileName);

% what is there?
quants = {VTKdata(:).datname};
[selection,ok] = listdlg('ListString',quants,'SelectionMode','single',...
       'PromptString','Select quantity to plot');
if(~ok)
    return
end

% extract data
quant = quants(selection);
vdata = VTKdata(strcmpi(quant,{VTKdata(:).datname})).datval;
qs = size(vdata);
qss = size(qs);

% get component to plot
if (qss(2)>3)
    qcomp = 0;
    while (qcomp==0)
        answer = inputdlg({['Select component to plot <= ' num2str(qs(1))]},...
            'Plot Component',1,{'1'});
        if(size(answer)==0)
            return
        end
        [val status] = str2num(answer{1});
        if status
            val = fix(val+.5);
            if (val>=1) && (val<=qs(1))
                qcomp = val;
            end
        end
    end
    % get 3D data array
    qdata = squeeze(vdata(qcomp,:,:,:));
else
    qdata = vdata;
end

% point coordinates
VTKfields = {VTKheader(:).parname};
origin = VTKheader(strcmpi('ORIGIN',VTKfields)).parval;
spacing = VTKheader(strcmpi('SPACING',VTKfields)).parval;
dims = size(qdata);
xloc = [origin(1):spacing(1):(origin(1)+(dims(1)-1)*spacing(1))];
yloc = [origin(2):spacing(2):(origin(2)+(dims(2)-1)*spacing(2))];
zloc = [origin(3):spacing(3):(origin(3)+(dims(3)-1)*spacing(3))];
[X,Y,Z] = ndgrid(xloc,yloc,zloc);

% vector list of points
npts = dims(1)*dims(2)*dims(3);
qv = reshape(qdata,npts,1);
xv = reshape(X,npts,1);
yv = reshape(Y,npts,1);
zv = reshape(Z,npts,1);

figure(1); scatter3(xv,yv,zv,3,qv,'filled')

