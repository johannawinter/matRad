function [Sigma_SISsq] = getSigmaSISsq(Identifier,basePath)

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
    

counter = 1;
i=1 ;
while(counter < 261)
  currentline = fgetl(fHandle);
   
        if(strfind(currentline,'energy'))
        energy(i)  = str2double(currentline(8:13));    
        focus = cell2mat(textscan(currentline,'energy %*f focus %f %f %f %f %f %f %f')); %at the HIT facility only beams with a FWHM >=6mm are used
        ix = find(focus>=6,1);
        FWHM(i) = focus(ix);
        i = i+1;
            
        end
        
    counter = counter +1;  
    
end

energy = energy';
FWHM = FWHM';
Sigma_SIS = FWHM / (2*sqrt(2*log(2)));
Sigma_SISsq = Sigma_SIS .* Sigma_SIS;

fclose(fHandle);


end
