function Pdata=js_getCarbon_Data

% extracts the first foci of each energy from the sis file
sigmaSISsq = tk_C_sigma_SISsq; 

for i = 1:255
    
    %Pdata(i).fileName = ['C_E' num2str(i) '_3MMRIFI_MS.ddd'];
    
    fid = fopen(['\\psf\Home\Documents\Heidelberg\TRiP98DATA\DDD\12C\RF3MM_NEU\12C_E' num2str(i) '_rifi3mm_ddd_new.ddd'],'r');
    
    if fid < 0
        display(['Could not open 12C_E' num2str(i) '_rifi3mm_ddd_new.ddd']);
    end
    % skip 7 lines:
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
% skip ende
    currentline = fgetl(fid);
    Pdata(i).energy = strread(currentline,'!energy %f');

%skip another 2 lines
fgetl(fid);
fgetl(fid);


depth = [];
ionization = [];
FWHM1 = [];
weight = [];
FWHM2 = [];

while ~feof(fid)
    depthn = fscanf(fid, '%f', 1);
    ionizationn = fscanf(fid, '%f',1);
    FWHM1n = fscanf(fid, '%f', 1);
    weightn = fscanf(fid, '%f', 1);
    FWHM2n = fscanf(fid, '%f', 1);
    fgetl(fid);
    
    depth = vertcat(depth, [depthn]);
    ionization = vertcat(ionization, [ionizationn]);
    FWHM1 =vertcat(FWHM1, [FWHM1n]);
    weight = vertcat(weight, [weightn]);
    FWHM2 = vertcat(FWHM2, [FWHM2n]);
end
    fclose(fid);    
   
    if i == 155
        st =2;
    end
    Pdata(i).depths = depth*10; 
    Pdata(i).Z   = ionization;
    %Pdata(i).DDD   = ionization;
    [val,idx]=max(ionization);
    Pdata(i).peakPos = Pdata(i).depths(idx)./10;
    %% prepare data for HIT
    h = sigmaSISsq( i,1);
    s11 = 1/(2*sqrt(2*log(2)))*abs(FWHM1);
    s22 = 1/(2*sqrt(2*log(2)))*abs(FWHM2);
    Pdata(i).sigma1 = sqrt( h + s11(:,1).^2);                                                                
    Pdata(i).sigma2 = sqrt( h + s22(:,1).^2);
    Pdata(i).weight = weight; 
    
%     Pdata(i).FWHM1 = FWHM1;
%     Pdata(i).FWHM2 = FWHM2;
%     Pdata(i).DDD   = ionization;
    % scaling for h5 file
%     [Pdata(i).absCalFac, ix] = max(Pdata(i).DDD); %find peak and peak-index ix
%     Pdata(i).peakPos = Pdata(i).depth(ix);
%     
%     Pdata(i).DDD = round(10000 * Pdata(i).DDD / Pdata(i).absCalFac); %DDD of each file Normed to peak-value * 10000
%     Pdata(i).absCalFac = Pdata(i).absCalFac ./ 10000;
    
end    