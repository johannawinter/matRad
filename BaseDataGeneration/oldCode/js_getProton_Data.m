function Pdata=js_getProton_Data
%% old function from Julian Streitz
for i = 1:255
    
    %Pdata(i).fileName = ['p_E' num2str(i) '_3MMRIFI_MS.ddd'];
    
    fid = fopen(['E:\TRiP98DATA_HIT-20131120\DDD\p\RF0MM\p_E' num2str(i) '_rifi0mm_ms.ddd'],'r');
    
    if fid < 0
        display(['Could not open p_E' num2str(i) '_rifi0mm_ms.ddd']);
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
    
    Pdata(i).depth = depth*10; 
    Pdata(i).Z   = ionization;
    [val,idx]=max(ionization);
    Pdata(i).peakPos = Pdata(i).depth(idx)./10;
    Pdata(i).sigma1 = FWHM1/(2*sqrt(2*log(2)));                                                               
    Pdata(i).sigma2 = FWHM2/(2*sqrt(2*log(2)));
    Pdata(i).weight = weight;   
     

    % scaling for h5 file
%     [Pdata(i).absCalFac, ix] = max(Pdata(i).DDD); %find peak and peak-index ix
%     Pdata(i).peakPos = Pdata(i).depth(ix);
%     
%     Pdata(i).DDD = round(10000 * Pdata(i).DDD / Pdata(i).absCalFac); %DDD of each file Normed to peak-value * 10000
%     Pdata(i).absCalFac = Pdata(i).absCalFac ./ 10000;
    
end    