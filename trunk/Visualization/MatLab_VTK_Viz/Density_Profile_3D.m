% MatLab Script: Density_Profile_3D.
%
% Date: 17 Feb 2011
% Author: John Nairn
%
% For any VTK file that include mass, this script will plot
% density profile along X, Y, or Z axis (by user input). The profile
% will be integrated within user-specified ranges for the other
% two axes.

% User select a file
[fileName,fldrPath,FilterIndex] = uigetfile('*.vtk','Select VTK File');
if(fileName==0)
    return
end

% read all data
[VTKheader,VTKdata] = ReadVTKFile(fldrPath,fileName,{'mass'});

% is it there?
qs=size(VTKdata);
if(qs(1)==0)
    h = msgbox('No mass data was found in that vtk file',...
        'Density Profile 3D','error','modal');
    return
end

% get mass data
qdata = VTKdata(1).datval;

% which direction ?
quants = {'X Axis','Y Axis','Z Axis'};
[axis,ok] = listdlg('ListString',quants,'SelectionMode','single',...
       'PromptString','Along which axis?');
if(~ok)
    return
end

% find x and y limits
VTKfields = {VTKheader(:).parname};
origin = VTKheader(strcmpi('ORIGIN',VTKfields)).parval;
spacing = VTKheader(strcmpi('SPACING',VTKfields)).parval;
dims = size(qdata);
xlim = [origin(1) (origin(1)+(dims(1)-1)*spacing(1))];
ylim = [origin(2) (origin(2)+(dims(2)-1)*spacing(2))];
zlim = [origin(3) (origin(3)+(dims(3)-1)*spacing(3))];
if (axis==1)
    prompt = {'Range within Y axis' 'Range within Z axis'};
    init = {num2str(ylim) num2str(zlim)};
elseif (axis==2)
    prompt = {'Range within X axis' 'Range within Z axis'};
    init = {num2str(xlim) num2str(zlim)};
else
    prompt = {'Range within X axis' 'Range within Y axis'};
    init = {num2str(xlim) num2str(ylim)};
end
hasLimits = 0;
while (hasLimits==0)
    answer = inputdlg(prompt,'Integration Range',1,init);
    if(size(answer)==0)
        return
    end
    [rng1 status1] = str2num(answer{1});
    [rng2 status2] = str2num(answer{2});
    if status1 && status2
        if size(rng1)==[1 2]
            if size(rng2)==[1 2]
                hasLimits=1;
            end
        end
    end
    if (hasLimits==0)
        h = msgbox('Each range must have two valid numbers',...
                        'Density Profile 3D','warn','modal');
        uiwait(h);
    end
end

% switch ranges if needed

% integrate area along axis

% plot the results

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

%figure(1); scatter3(xv,yv,zv,3,qv,'filled')

