function [C_sigma_SISsq] = tk_C_sigma_SISsq

fHandle = fopen('\\psf\Home\Documents\Heidelberg\TRiP98DATA\SIS\12C_1.6.2008.sis','r');


    if fHandle < 0
        display('Could not open SIS file');
    end
    
% fgetl(fHandle);
% fgetl(fHandle);   
% fgetl(fHandle);
% fgetl(fHandle);
% fgetl(fHandle);
% fgetl(fHandle);
% 
% currentline = fgetl(fHandle);
% energy = str2double(currentline(8:13));


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
C_sigma_SIS = FWHM / (2*sqrt(2*log(2)));
C_sigma_SISsq = C_sigma_SIS .* C_sigma_SIS;
energy_C_sigma_SISsq = [energy C_sigma_SISsq];


fclose(fHandle);


end
