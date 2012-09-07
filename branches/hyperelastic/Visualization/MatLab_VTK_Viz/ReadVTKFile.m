function [ VTKheader,VTKdata ] = ReadVTKFile( fldrPath,fileName,wantData )
% fldrPath is path to folder containing file
%    (or empty to use current folder)
% fileName for vtk file (concatonated to flrdPath to get full path)
% wantData is cell array of data types to return as named in vtk file
%    (e.g., strain, stress, totalstrain, etc.) (optional)

% create wantData if needed
if(nargin<3)
    wantData = {};
    maxData = 10000;
else
    csize = size(wantData);
    maxData = csize(2);
end

% Name of fields in header that are read in this function
VTKfields={'step';'DIMENSIONS';'ORIGIN';'SPACING'};

% Initiate dataset paremeter structure file
VTKheader=struct('parname',VTKfields, 'parval',[]);

% Initiate dataset 
VTKdata = struct('datname',{},'datval',[]); 

% Initiate dataset number and the header line number
dsn = 0;

% Open file and read first line
fid = fopen([fldrPath fileName],'r');
S = fgetl(fid);

% Loop over entire file
while ~isnumeric(S)
    parname = sscanf(S,'%s',1);
    
    %filter out important parameters of the VTK file header and data
    switch parname
        case {'DIMENSIONS','ORIGIN','SPACING'}
            VTKheader(strcmpi(parname,VTKfields)).parval = ...
                sscanf(S,'%*s %f %f %f');
        case 'POINT_DATA'
            %do nothing
        case 'TENSORS'
            vname = sscanf(S,'%*s %s',1);
            sizeVTK = [9; VTKheader(strcmpi('DIMENSIONS',VTKfields)).parval]';
            vdata=reshape(fscanf(fid,'%f',prod(sizeVTK)),sizeVTK);
            if(and(nargin>2,strcmpi(vname,wantData)==0))
                vdata = [];
            else
                dsn = dsn +1; %advance the dataset number
                VTKdata(dsn).datname = vname;
                VTKdata(dsn).datval = vdata;
                if(dsn==maxData)
                    break;
                    
                end
            end
        case 'SCALARS'
            vname = sscanf(S,'%*s %s',1);
            sizeVTK = [VTKheader(strcmpi('DIMENSIONS',VTKfields)).parval]';
            %pass the 'LOOKUP_TABLE' line
            S = fgetl(fid);
            vdata=reshape(fscanf(fid,'%f',prod(sizeVTK)),sizeVTK);
            if(and(nargin>2,strcmpi(vname,wantData)==0))
                vdata = [];
            else
                dsn = dsn +1; %advance the dataset number
                VTKdata(dsn).datname = vname;
                VTKdata(dsn).datval = vdata;
                if(dsn==maxData)
                    break;
                end
            end
        case 'VECTORS'
            vname = sscanf(S,'%*s %s',1);
            sizeVTK = [3; VTKheader(strcmpi('DIMENSIONS',VTKfields)).parval]'
            vdata=reshape(fscanf(fid,'%f',prod(sizeVTK)),sizeVTK);
            if(and(nargin>2,strcmpi(vname,wantData)==0))
                vdata = [];
            else
                dsn = dsn +1; %advance the dataset number
                VTKdata(dsn).datname = vname;
                VTKdata(dsn).datval = vdata;
                if(dsn==maxData)
                    break;
                end
            end
        case 'LOOKUP_TABLE'
            %do nothing
        case 'ASCII'
            %do nothing
        case 'DATASET'
            %do nothing
        case '#'
            %do nothing
        otherwise
            k = strfind(parname,'step');
            if (k==1)
                VTKheader(strcmpi('step',VTKfields)).parval = ...
                    sscanf(S,'step:%f time: %f');
            end
    end
    
    %scan next line  
    S = fgetl(fid);
end

fclose(fid);

end

