function ddd=matRad_getDDDfromTxt(Identifier,basePath,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_getDDDfromTxt script
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
% This script can be used to parse depth dose distribution from txt files.
% The identifier is important to chose the correct relative path.
% if double lateral sigmas are available in the ddd file then they will be
% used otherwise constant sigma or the highland formula

switch Identifier
    case 'C'
        relPath = [filesep 'DDD' filesep '12C' filesep 'RF3MM_NEU'];
        %relPath = [filesep 'DDD\12C\RF3MM'];
    case {'p','H'}
        relPath = [filesep 'DDD' filesep' 'p' filesep 'RF0MM'];
    case 'O'
        relPath = [filesep 'DDD' filesep '16O' filesep 'RF3MM'];
    otherwise
        error('unkown particle type')
end

% try to read inital widths of the beam
try
    sigmaSISsq = getSigmaSISsq(Identifier,basePath);
catch
    warning('Couldnt read sis files containing the inital widths');
    sigmaSISsq = 0;
end

Files = dir([basePath relPath]);

for k = length(Files):-1:1
    % remove hidden folders starting with . and .. or .DSStore
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

    % allocate some variables
    depth = [];
    ionization = [];
    FWHM1 = [];
    weight = [];
    FWHM2 = [];
    
    % parse line
    while ~feof(fid)
        
        LineData = str2num(fgetl(fid));
        depth = vertcat(depth, [LineData(1,1)]);
        ionization = vertcat(ionization, [LineData(1,2)]);
        
        if length(LineData)==5
            FWHM1 =vertcat(FWHM1, [LineData(1,3)]);
            weight = vertcat(weight, [LineData(1,4)]);
            FWHM2 = vertcat(FWHM2, [LineData(1,5)]);
        end
    end
    fclose(fid);    

    ddd(i).depths = depth*10; 
    ddd(i).Z   = ionization;
    [val,idx]=max(ionization);
    ddd(i).peakPos = ddd(i).depths(idx);
    
    %% calculate sigmas and weights used to model lateral profile
    if ~isempty(FWHM1)
        s11 = 1/(2*sqrt(2*log(2)))*abs(FWHM1);
        s22 = 1/(2*sqrt(2*log(2)))*abs(FWHM2);
        ddd(i).sigma1 = sqrt( sigmaSISsq( i,1) + s11(:,1).^2);                                                                
        ddd(i).sigma2 = sqrt( sigmaSISsq( i,1) + s22(:,1).^2);
        ddd(i).weight = weight; 
    else
        %% TODO add analytical calculation of sigma
        ddd(i).sigma = ones(size(ddd(i).depths,1),1)*sqrt(sigmaSISsq(i,1));
    end

    
    
end   

% sort content according to enertgy;
[~,IdxEnergy] = sort([ddd.energy]);
ddd=ddd(IdxEnergy);

%% plot random ddd 
if visBool
    
    switch Identifier
        case 'C'
             load carbonBaseData
        case {'p','H'}
            load protonBaseData       
    end
    
    vIdx = round(linspace(10,length(baseData)-10,4));
    figure,
    
    for i = 1:length(vIdx)
        Energy = baseData(vIdx(i)).energy;
        [~,IdxDDD]=min(abs([ddd.energy]-Energy));
        subplot(2,2,i),plot(ddd(IdxDDD).depths,ddd(IdxDDD).Z,'r'),hold on
                       plot(baseData(vIdx(i)).depths,baseData(vIdx(i)).Z,'b');
                       legend({'new ddd - HIT','existing baseData'});
                       xlabel(' depth in mm '), ylabel('Z in Gy'), title(['Energy: ' num2str(Energy)])
    end
    
    
    
end