function mx = matRad_readSparseMx(filename)

% open file for writing
h = fopen(filename,'rb');

% number of rows
numOfRows = fread(h,1,'integer*4');

% number of columns
numOfColumns = fread(h,1,'integer*4');

% number of nonzeros
numOfNonzeros = fread(h,1,'integer*4');

% write rows
i = fread(h,numOfNonzeros,'integer*4');

% write columns
j = fread(h,numOfNonzeros,'integer*4');

% write values
s = fread(h,numOfNonzeros,'double');

% close file handle
fclose(h);

% find row and coumn indices including values
i(1000)
j(1000)
s(1000)

% generate sparse matrix
mx = sparse(i,j,s,numOfRows,numOfColumns,numOfNonzeros);