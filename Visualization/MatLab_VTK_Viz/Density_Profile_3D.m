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
if (rng1(1)>rng1(2))
    rng1 = [rng1(2) rng1(1)];
end
if (rng2(1)>rng2(2))
    rng2 = [rng2(2) rng2(1)];
end

% integrate area along axis
if(axis==1)
    lim1 = [max(1,fix((rng1(1)-origin(2))/spacing(2))):...
                 min(dims(2),fix((rng1(2)-origin(2))/spacing(2)))];
    lim2 = [max(1,fix((rng2(1)-origin(3))/spacing(3))):...
                 min(dims(3),fix((rng2(2)-origin(3))/spacing(3)))];
else
    lim1 = [max(1,fix((rng1(1)-origin(1))/spacing(1))):...
                 min(dims(1),fix((rng1(2)-origin(1))/spacing(1)))];
    if(axis==3)
        lim2 = [max(1,fix((rng2(1)-origin(2))/spacing(2))):...
                 min(dims(2),fix((rng2(2)-origin(2))/spacing(2)))];
    else
        lim2 = [max(1,fix((rng2(1)-origin(3))/spacing(3))):...
                 min(dims(3),fix((rng2(2)-origin(3))/spacing(3)))];
    end
end

tdat = [];
ddat = [];
if(axis==1)
    for i = 1:dims(1)
        tdat(i) = origin(1)+(i-1)*spacing(1);
        dplane = squeeze(qdata(i,lim1,lim2));
        ddat(i) = sum(sum(dplane));
    end
elseif(axis==2)
    for i = 1:dims(2)
        tdat(i) = origin(2)+(i-1)*spacing(2);
        dplane = squeeze(qdata(lim1,i,lim2));
        ddat(i) = sum(sum(dplane));
    end
else
    for i = 1:dims(3)
        tdat(i) = origin(3)+(i-1)*spacing(3);
        dplane = squeeze(qdata(lim1,lim2,i));
        ddat(i) = sum(sum(dplane));
    end
end

% scale to denity in g/cm^3
s1 = size(lim1);
s2 = size(lim2);
cellVolume = (s1(2)-1)*(s2(2)-1)*spacing(1)*spacing(2)*spacing(3)/1000.;
ddat = ddat/cellVolume;

% print average density
ds = size(ddat);
avgDensity = sum(ddat)/ds(2);

% plot the results
figure(1);
p = plot(tdat,ddat);
xlabel('Position (mm)')
ylabel('Density (g/cm^3)')
title('MPM Density Profile')
text(tdat(fix(ds(2)/2)),.05,['avg = ' num2str(avgDensity)])

