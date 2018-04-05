function nrrd = matRad_readRegisteredNrrd(filename,ct)

% copy some variables
nrrd.x          = ct.x;
nrrd.y          = ct.y;
nrrd.z          = ct.z;
nrrd.resolution = ct.resolution;

% read cube
[nrrd_orig, meta] = nrrdread(filename);

% get some meta data from nrrd cube
[xOrigin, yOrigin, zOrigin] = strread(meta.spaceorigin,'(%f,%f,%f)');
[xDim,    yDim,    zDim]    = strread(meta.sizes,'%f %f %f');
[xRes,    yRes,    zRes]    = strread(meta.spacedirections,'(%f,0,0) (0,%f,0) (0,0,%f)');

% set up grid vectors for original nrrd cube
%xCoords = xOrigin + xRes * [0:xDim-1];
%yCoords = yOrigin + yRes * [0:yDim-1];
xCoords = ct.dicomInfo.ImagePositionPatient(1) + xRes * [0:xDim-1];
yCoords = ct.dicomInfo.ImagePositionPatient(2) + yRes * [0:yDim-1];
zCoords = zOrigin + zRes * [0:zDim-1];

% cubic interpolation (extrapolation yields 0)
nrrd.cube = interp3(xCoords,yCoords',zCoords,nrrd_orig, ...
                    nrrd.x,nrrd.y',nrrd.z,'cubic',0);

for i = 1:size(ct.cube,3)
    subplot(1,2,1)
    imagesc(ct.cube(:,:,i))
    colormap gray
    colorbar
    title(['original ct - slice ' num2str(i)])
    subplot(1,2,2)
    imagesc(nrrd.cube(:,:,i))
    colormap gray
    colorbar
    title(['interpolated from nrrd - slice ' num2str(i)])
    drawnow
    pause(.05)
end

end


function [X, meta] = nrrdread(filename)
%NRRDREAD  Import NRRD imagery and metadata.
%   [X, META] = NRRDREAD(FILENAME) reads the image volume and associated
%   metadata from the NRRD-format file specified by FILENAME.
%
%   Current limitations/caveats:
%   * "Block" datatype is not supported.
%   * Only tested with "gzip" and "raw" file encodings.
%   * Very limited testing on actual files.
%   * I only spent a couple minutes reading the NRRD spec.
%
%   See the format specification online:
%   http://teem.sourceforge.net/nrrd/format.html

% Copyright 2012 The MathWorks, Inc.


% Open file.
fid = fopen(filename, 'rb');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% Magic line.
theLine = fgetl(fid);
assert(numel(theLine) >= 4, 'Bad signature in file.')
assert(isequal(theLine(1:4), 'NRRD'), 'Bad signature in file.')

% The general format of a NRRD file (with attached header) is:
% 
%     NRRD000X
%     <field>: <desc>
%     <field>: <desc>
%     # <comment>
%     ...
%     <field>: <desc>
%     <key>:=<value>
%     <key>:=<value>
%     <key>:=<value>
%     # <comment>
% 
%     <data><data><data><data><data><data>...

meta = struct([]);

% Parse the file a line at a time.
while (true)

  theLine = fgetl(fid);
  
  if (isempty(theLine) || feof(fid))
    % End of the header.
    break;
  end
  
  if (isequal(theLine(1), '#'))
      % Comment line.
      continue;
  end
  
  % "fieldname:= value" or "fieldname: value" or "fieldname:value"
  parsedLine = regexp(theLine, ':=?\s*', 'split','once');
  
  assert(numel(parsedLine) == 2, 'Parsing error')
  
  field = lower(parsedLine{1});
  value = parsedLine{2};
  
  field(isspace(field)) = '';
  meta(1).(field) = value;
  
end

datatype = getDatatype(meta.type);

% Get the size of the data.
assert(isfield(meta, 'sizes') && ...
       isfield(meta, 'dimension') && ...
       isfield(meta, 'encoding') && ...
       isfield(meta, 'endian'), ...
       'Missing required metadata fields.')

dims = sscanf(meta.sizes, '%d');
ndims = sscanf(meta.dimension, '%d');
assert(numel(dims) == ndims);

data = readData(fid, meta, datatype);
data = adjustEndian(data, meta);

% Reshape and get into MATLAB's order.
X = reshape(data, dims');
X = permute(X, [2 1 3]);

end

function datatype = getDatatype(metaType)

% Determine the datatype
switch (metaType)
 case {'signed char', 'int8', 'int8_t'}
  datatype = 'int8';
  
 case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
  datatype = 'uint8';

 case {'short', 'short int', 'signed short', 'signed short int', ...
       'int16', 'int16_t'}
  datatype = 'int16';
  
 case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
       'uint16_t'}
  datatype = 'uint16';
  
 case {'int', 'signed int', 'int32', 'int32_t'}
  datatype = 'int32';
  
 case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
  datatype = 'uint32';
  
 case {'longlong', 'long long', 'long long int', 'signed long long', ...
       'signed long long int', 'int64', 'int64_t'}
  datatype = 'int64';
  
 case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
       'uint64', 'uint64_t'}
  datatype = 'uint64';
  
 case {'float'}
  datatype = 'single';
  
 case {'double'}
  datatype = 'double';
  
 otherwise
  assert(false, 'Unknown datatype')
end

end

function data = readData(fidIn, meta, datatype)

switch (meta.encoding)
 case {'raw'}
  
  data = fread(fidIn, inf, [datatype '=>' datatype]);
  
 case {'gzip', 'gz'}

  tmpBase = tempname();
  tmpFile = [tmpBase '.gz'];
  fidTmp = fopen(tmpFile, 'wb');
  assert(fidTmp > 3, 'Could not open temporary file for GZIP decompression')
  
  tmp = fread(fidIn, inf, 'uint8=>uint8');
  fwrite(fidTmp, tmp, 'uint8');
  fclose(fidTmp);
  
  gunzip(tmpFile)
  
  fidTmp = fopen(tmpBase, 'rb');
  cleaner = onCleanup(@() fclose(fidTmp));
  
  meta.encoding = 'raw';
  data = readData(fidTmp, meta, datatype);
  
 case {'txt', 'text', 'ascii'}
  
  data = fscanf(fidIn, '%f');
  data = cast(data, datatype);
  
 otherwise
  assert(false, 'Unsupported encoding')
end

end

function data = adjustEndian(data, meta)

[~,~,endian] = computer();

needToSwap = (isequal(endian, 'B') && isequal(lower(meta.endian), 'little')) || ...
             (isequal(endian, 'L') && isequal(lower(meta.endian), 'big'));
         
if (needToSwap)
    data = swapbytes(data);
end

end

