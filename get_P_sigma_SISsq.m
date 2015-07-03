function [P_sigma_SISsq] = get_P_sigma_SISsq

fHandle = fopen('E:\TRiP98DATA_HIT-20131120\SIS\1H_1.1.2009.sis','r');


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
        FWHM(i) = str2double(currentline(20:24)); 
        i = i+1;
            
        end
        
    counter = counter +1;
    
    
    
end

energy = energy';
FWHM = FWHM';
P_sigma_SIS = FWHM / (2*sqrt(2*log(2)));
P_sigma_SISsq = P_sigma_SIS .* P_sigma_SIS;
energy_P_sigma_SISsq = [energy P_sigma_SISsq];














fclose(fHandle);




end
