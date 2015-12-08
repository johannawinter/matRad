function matRad_writeSparseMx(mx,filename)

% open file for writing
h = fopen(filename,'wb');

% number of rows
fwrite(h,size(mx,1),'integer*4');

% number of columns
fwrite(h,size(mx,2),'integer*4');

% number of nonzeros
fwrite(h,nnz(mx),'integer*4');

% find row and coumn indices including values
[i,j,s] = find(mx);

i(1000)
j(1000)
s(1000)

% write rows
fwrite(h,i,'integer*4');

% write columns
fwrite(h,j,'integer*4');

% write values
fwrite(h,s,'double');

% close file handle
fclose(h);