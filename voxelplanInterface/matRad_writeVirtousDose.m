function matRad_writeVirtousDose(resultGUI,ct,filename)

% interpolate dose
xi = ([1:ct.virtuosHeader.dimy] -.5) * ct.virtuosHeader.pixel_size;
yi = ([1:ct.virtuosHeader.dimx]' -.5) * ct.virtuosHeader.pixel_size;
zi = ct.virtuosHeader.z_table.position;

d = interp3(ct.x,ct.y,ct.z,resultGUI.physicalDose,xi,yi,zi);

% rotate
for i = 1:size(d,3)
    d(:,:,i) = rot90(flipud(d(:,:,i)),-1);
end
    
% write data
h = fopen([filename '.dos'],'w');
fwrite(h,d(:),'single');
fclose(h);