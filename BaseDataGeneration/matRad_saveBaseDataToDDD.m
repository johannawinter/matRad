load('..\protonBaseData')

metaInfo.filetype    = 'ddd';
metaInfo.fileversion = '20151610';
metaInfo.filedate    = 'Tue Oct 16 12:00:00 2015';
metaInfo.projectile  =  '1H';
metaInfo.material    =  'H20';
metaInfo.composition =  'H20';
metaInfo.generatedBy =  'Hans-Peter Wieser';
metaInfo.density     =  '1';


fNames = fieldnames(metaInfo);
for i = 1:length(baseData)
    FileName = ['p_E' num2str(i) '_rifi0mm_ms.ddd'];
    fileID = fopen(FileName,'w');
    for j = 1:length(fNames)
        fprintf(fileID,['!' fNames{j,1} '    ' metaInfo.(fNames{j,1}) '\n']);
    end
    
    fprintf(fileID,['!energy ' num2str(baseData(i).energy) '\n']);
    fprintf(fileID,'#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] sigma[mm]\n');
    fprintf(fileID,'!ddd\n');
    
    A = [baseData(i).depths./10, baseData(i).Z,baseData(i).sigma]';
    fprintf(fileID,'%5.4f %10.7f %5.4f\r\n',A);
    fclose(fileID);
end
clear metaInfo


%%
load('..\carbonBaseData')

metaInfo.filetype    = 'ddd';
metaInfo.fileversion = '20151610';
metaInfo.filedate    = 'Tue Oct 16 12:00:00 2015';
metaInfo.projectile  =  '1H';
metaInfo.material    =  'H20';
metaInfo.composition =  'H20';
metaInfo.generatedBy =  'Hans-Peter Wieser';
metaInfo.density     =  '1';


fNames = fieldnames(metaInfo);
for i = 1:length(baseData)
    FileName = ['12C_E' num2str(i) '_rifi3mm_ms.ddd'];
    fileID = fopen(FileName,'w');
    for j = 1:length(fNames)
        fprintf(fileID,['!' fNames{j,1} '    ' metaInfo.(fNames{j,1}) '\n']);
    end
    
    fprintf(fileID,['!energy ' num2str(baseData(i).energy) '\n']);
    fprintf(fileID,'#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] sigma[mm]\n');
    fprintf(fileID,'!ddd\n');
    
    A = [baseData(i).depths./10, baseData(i).Z,baseData(i).sigma]';
    fprintf(fileID,'%5.4f %10.7f %5.4f\r\n',A);
    fclose(fileID);
end



