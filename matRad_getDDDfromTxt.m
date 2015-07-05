function ddd=matRad_getDDDfromTxt(Identifier,basePath)

% extracts the first foci of each energy from the sis file
% if other foci shold be used, adapt the function

switch Identifier
    case 'C'
        relPath = [filesep 'DDD\12C\RF3MM_NEU'];
    case {'p','H'}
        relPath = [filesep 'DDD\p\RF0MM'];
    case 'O'
        relPath = [filesep 'DDD\16O\RF3MM'];
    otherwise
        error('unkown particle type')
end

sigmaSISsq = getSigmaSISsq(Identifier,basePath);
Files = dir([basePath relPath]);

for k = length(Files):-1:1
    % remove hidden folders starting with .
    fname = Files(k).name;
    if fname(1) == '.'
        Files(k) = [ ];
    end
end

for i = 1:length(Files)
    
    fid = fopen([basePath relPath filesep Files(i).name],'r');
    if fid < 0
        display(['Could not open ' path filesep Files(i).name]);
    end
    % skip 7 lines:
    for j= 1:7
        fgetl(fid);
    end

    % read energy
    currentline = fgetl(fid);
    ddd(i).energy = strread(currentline,'!energy %f');

    %skip another 2 lines
    fgetl(fid);
    fgetl(fid);


    depth = [];
    ionization = [];
    FWHM1 = [];
    weight = [];
    FWHM2 = [];
    
    % parse line
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

    ddd(i).depths = depth*10; 
    ddd(i).Z   = ionization;
    [val,idx]=max(ionization);
    ddd(i).peakPos = ddd(i).depths(idx)./10;
    
    %% calculate sigmas and weights used to model lateral profile
    s11 = 1/(2*sqrt(2*log(2)))*abs(FWHM1);
    s22 = 1/(2*sqrt(2*log(2)))*abs(FWHM2);
    ddd(i).sigma1 = sqrt( sigmaSISsq( i,1) + s11(:,1).^2);                                                                
    ddd(i).sigma2 = sqrt( sigmaSISsq( i,1) + s22(:,1).^2);
    ddd(i).weight = weight; 
    
    
end   

% sort content according to enertgy;
[~,IdxEnergy] = sort([ddd.energy]);
ddd=ddd(IdxEnergy);


