function matRad_writeVirtousDose(resultGUI,ct,patFolder,filename)

fprintf('write *.dos and *.hed file\n');

%% write *.dos file
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
h = fopen([patFolder filesep filename '.dos'],'w');
fwrite(h,d(:),'single');
fclose(h);

%% write corresponding header
h = fopen([patFolder filesep filename '.hed'],'w');
fprintf(h,'version 1.4\n');
fprintf(h,'modality ct\n');
fprintf(h,'created_by \n');
fprintf(h,'creation_info \n');
fprintf(h,'primary_view transversal\n');
fprintf(h,'data_type float\n');
fprintf(h,'num_bytes 4\n');
fprintf(h,'byte_order linux\n');
fprintf(h,['patient_name ' filename '\n']);
fprintf(h,'slice_dimension %d\n',ct.virtuosHeader.dimx);
fprintf(h,'pixel_size %f\n',ct.virtuosHeader.pixel_size);
fprintf(h,'slice_distance %f\n',ct.virtuosHeader.slice_distance);
fprintf(h,'slice_number %d\n',ct.virtuosHeader.dimz);
fprintf(h,'xoffset 0\n');
fprintf(h,'dimx %d\n',ct.virtuosHeader.dimx);
fprintf(h,'yoffset 0\n');
fprintf(h,'dimy %d\n',ct.virtuosHeader.dimx);
fprintf(h,'zoffset 0\n');
fprintf(h,'dimz %d\n',ct.virtuosHeader.dimz);
fprintf(h,'z_table yes\n');
fprintf(h,'slice_no  position  thickness  gantry_tilt\n');
for i = 1:ct.virtuosHeader.dimz
    fprintf(h,'%d %f %f %f\n',ct.virtuosHeader.z_table.slice_no(i), ...
                              ct.virtuosHeader.z_table.position(i), ...
                              ct.virtuosHeader.z_table.thickness(i), ...
                              ct.virtuosHeader.z_table.gantry_tilt(i));
end
fclose(h);