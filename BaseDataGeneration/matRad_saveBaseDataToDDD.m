clc, clear 
close
load('..\protons_HIT')

metaInfo.filetype    = 'ddd';
metaInfo.fileversion = '20151610';
metaInfo.filedate    = 'Tue Nov 11 12:00:00 2015';
metaInfo.projectile  =  '1H';
metaInfo.material    =  'H20';
metaInfo.composition =  'H20';
metaInfo.generatedBy =  'Hans-Peter Wieser';
metaInfo.density     =  '1';
metaInfo.description =  'ddd from Trip98, lateral info from kata.p - inital beam width is not considered';

Factor = 2*sqrt(2*log(2));
[Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS('p','E:\TRiP98DATA_HIT-20131120',1);

fNames = fieldnames(metaInfo);
for i = 1:length(machine.data)
    FileName = ['p_E' num2str(i) '_rifi0mm_ms.ddd'];
    fileID = fopen(FileName,'w');
    for j = 1:length(fNames)
        fprintf(fileID,['!' fNames{j,1} '    ' metaInfo.(fNames{j,1}) '\n']);
    end
    
    fprintf(fileID,['!energy ' num2str(machine.data(i).energy) '\n']);
    fprintf(fileID,'#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[mm] weight[a.u.] FWHM2[mm]\n');
    fprintf(fileID,'!ddd\n');
    
    A = [machine.data(i).depths./10, machine.data(i).Z,sqrt(machine.data(i).sigma1(:).^2 - Sigma_SIS(i).^2).*Factor,...
        machine.data(i).weight,   sqrt(machine.data(i).sigma2(:).^2 - Sigma_SIS(i)^2).*Factor]';
    fprintf(fileID,'%5.4f %10.7f %5.4f %5.4f %5.4f \r\n',A);
    fclose(fileID);
end
clear metaInfo


%%
load('..\carbon_HIT')
[Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS('C','E:\TRiP98DATA_HIT-20131120',1);
metaInfo.filetype    = 'ddd';
metaInfo.fileversion = '20151610';
metaInfo.filedate    = 'Tue Nov 25 12:00:00 2015';
metaInfo.projectile  =  '12C';
metaInfo.material    =  'H20';
metaInfo.composition =  'H20';
metaInfo.generatedBy =  'Hans-Peter Wieser';
metaInfo.density     =  '1';
metaInfo.description =  'ddd from Trip98, lateral info from kata.p - inital beam width is not considered';


fNames = fieldnames(metaInfo);
for i = 1:length(machine.data)
    FileName = ['12C_E' num2str(i) '_rifi3mm_ms.ddd'];
    fileID = fopen(FileName,'w');
    for j = 1:length(fNames)
        fprintf(fileID,['!' fNames{j,1} '    ' metaInfo.(fNames{j,1}) '\n']);
    end
    
    fprintf(fileID,['!energy ' num2str(machine.data(i).energy) '\n']);
   fprintf(fileID,'#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[mm] weight[a.u.] FWHM2[mm]\n');
    fprintf(fileID,'!ddd\n');
    
 A = [machine.data(i).depths./10, machine.data(i).Z,sqrt(machine.data(i).sigma1(:).^2 - Sigma_SIS(i).^2).*Factor,...
        machine.data(i).weight,   sqrt(machine.data(i).sigma2(:).^2 - Sigma_SIS(i)^2).*Factor]';
 fprintf(fileID,'%5.4f %10.7f %5.4f %5.4f %5.4f \r\n',A); 
    fclose(fileID);
end



