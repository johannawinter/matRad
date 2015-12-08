function [ct,cst] = matRad_importVirtuosDataSet(patientID,writeDataBool)

if nargin < 2
    writeDataBool = 0;
end

%% header
fprintf('read hed file...\n');
headerHandle = fopen([patientID '000.hed'],'r');

if headerHandle < 0
    fprintf('Failed to open .hed file.\n')
    return;
end

header.z_table.slice_no    = [];
header.z_table.position    = [];
header.z_table.thickness   = [];
header.z_table.gantry_tilt = [];

while ~feof(headerHandle)
    
    currentLine = fgetl(headerHandle);

     if strfind(currentLine,'pixel_size')
        header.pixel_size      = str2double(currentLine(12:end));
     elseif strfind(currentLine,'slice_distance')
        header.slice_distance  = str2double(currentLine(16:end));
     elseif strfind(currentLine,'slice_number')
        header.slice_number    = str2double(currentLine(14:end));
     elseif strfind(currentLine,'dimx')
        header.dimx            = str2double(currentLine(6:end));
     elseif strfind(currentLine,'dimy')
        header.dimy            = str2double(currentLine(6:end));
     elseif strfind(currentLine,'dimz')
        header.dimz            = str2double(currentLine(6:end));
     elseif strfind(currentLine,'z_table')
        header.z_table.bool   = currentLine(9:end);
        % skip one line
        currentLine = fgetl(headerHandle);
        tmp = cell2mat(textscan(headerHandle,'%f %f %f %f',header.dimz));
        header.z_table.slice_no    = tmp(:,1);
        header.z_table.position    = tmp(:,2);
        header.z_table.thickness   = tmp(:,3);
     end   
     
end

fclose(headerHandle);

%% ct
% read data from file
fprintf('read ctx file...\n');
ctHandle = fopen([patientID '000.ctx'],'r');
if ctHandle < 0
    error('Could not load ct.\n')
end
ctVector = fread(ctHandle,inf, 'short');
fclose(ctHandle);

if header.dimx*header.dimy*header.dimz ~= size(ctVector)
    error('Dimensions of dose cube given in .hed file and in .cti/.ctx file differ.')
else
    origCt = double(reshape(ctVector,header.dimx,header.dimy,header.dimz));
    for i = 1:size(origCt,3)
        origCt(:,:,i) = flipud(rot90(squeeze(origCt(:,:,i))));
    end
end  

% interpolation
resolution.x = 2.5;
resolution.y = 2.5;
resolution.z = 2.5;

info.ImagePositionPatient = [.5*header.pixel_size*ones(1,header.slice_number); ...
                             .5*header.pixel_size*ones(1,header.slice_number); ...
                             header.z_table.position'];
                         
info.Columns      = header.dimx;
info.Rows         = header.dimy;
info.PixelSpacing = [header.pixel_size header.pixel_size header.slice_distance];

% interpolate ct
fprintf('interpolate ct...\n');
ct = matRad_interpCtCube(origCt, info, resolution);

% conversion from HU to e- density
fprintf('conversion of HU to waterEqD...\n');
load hlutDefault.mat; % load LUT

% Manual adjustments if ct data is corrupt. If some values are out of range
% of the LUT, then these values are adjusted.
if sum(ct.cube > max(hlut(:,1)) | ct.cube < min(hlut(:,1)))>0
    warning('projecting out of range HU values');
    ct.cube(ct.cube<min(hlut(:,1))) = min(hlut(:,1));
    ct.cube(ct.cube>max(hlut(:,1))) = max(hlut(:,1));
end

% interpolate HU to relative electron density based on lookup table
ct.cube = interp1(hlut(:,1),hlut(:,2),double(ct.cube));

% save hlut
ct.hlut = hlut;

%% parse vdx
fprintf('parsing vdx file...\n');
% open file for reading
vdxHandle = fopen([patientID '000.vdx'],'r');
if vdxHandle < 0
    error('Could not open .vdx file. Aborting...\n')
end
% go through file line by line
structCounter = 0;
while ~feof(vdxHandle)    
    currentLine = fgetl(vdxHandle);
    % if key word voi is found -> initialze new voi
    if strncmp(currentLine,'voi ',4)
        structCounter = structCounter + 1;
        structures(structCounter).structName = cell2mat(strread(currentLine,'voi %s'));
        itemCounter = 0;
    end
    
    if strncmp(currentLine,'number_of_points ',17)
        itemCounter = itemCounter + 1;
        numOfPoints = strread(currentLine,'number_of_points %f');
        tmp = cell2mat(textscan(vdxHandle,'%f %f %f %f %f %f',numOfPoints));
        structures(structCounter).item(itemCounter).points = tmp(:,1:3);
    end
end

% prepare structures to use matRad dicom import functions
ct.dicomInfo.SlicePositions = header.z_table.position';
ct.dicomInfo.SliceThickness = header.z_table.thickness';

for i = 1:size(structures,2)
        fprintf('creating cube for %s volume...\n',structures(i).structName);
        structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
end

%% create cst struct
fprintf('create cst struct...\n');
cst = matRad_createCst(structures);

% add virtuos meta info about orig ct to ct sturct
ct.virtuosHeader = header;

if writeDataBool
    % write file
    fprintf('write file to disk...\n');
    save(patientID,'ct','cst');
end
