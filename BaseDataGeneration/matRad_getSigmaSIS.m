function [Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,basePath,FocusIdx)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getSigmaSISsq script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script can be used to parse the first foci of a sis file. Please
% note that it returns the sigma squared

switch Identifier
    case 'C'
        FullFilePath = [basePath filesep 'SIS' filesep '12C_1.6.2008.sis'];
    case {'p','H'}
        FullFilePath = [basePath filesep 'SIS' filesep '1H_1.1.2009.sis'];
    case 'O'
        FullFilePath = [basePath filesep 'SIS' filesep '16O_20.12.2010'];
    otherwise
        error('unkown particle type')
end


fHandle = fopen(FullFilePath,'r');

if fHandle < 0
    display('Could not open SIS file');
end
    
Cnt=1;
currentline = fgetl(fHandle);
while ischar(currentline)

    if(strfind(currentline,'energy'))
        vEnergySIS(Cnt,1)  = str2double(currentline(8:13));    
        vFoci = cell2mat(textscan(currentline,'energy %*f focus %f %f %f %f %f %f %f')); 
        FWHM(Cnt) = vFoci(FocusIdx);
        Cnt = Cnt+1;
    end
    currentline = fgetl(fHandle);       
end

Sigma_SIS = (FWHM./(2*sqrt(2*log(2))))';

fclose(fHandle);

end
