function ddd=matRad_getDDDfromTxt(Identifier,basePath,FocusIdx,visBool)
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

% try to read inital beam width which is energy and and machine specific
try
    sigmaSISsq = matRad_getSigmaSISsq(Identifier,basePath,FocusIdx);
catch
    warning('Could not read sis files containing the inital widths');
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
    % for basedata having single lateral gaussian
    FWHM = [];
    % for basedata having double lateral gaussian
    FWHM1 = [];
    weight = [];
    FWHM2 = [];
    
    % parse line
    currentLine = fgetl(fid);
    
    while ischar(currentLine)
        currentLine = str2num(currentLine);
        depth       = vertcat(depth,currentLine(1));
        ionization  = vertcat(ionization,currentLine(2));
        
        if length(currentLine)==5
            FWHM1  = vertcat(FWHM1,currentLine(3));
            weight = vertcat(weight,currentLine(1,4));
            FWHM2  = vertcat(FWHM2,currentLine(1,5));
        elseif length(currentLine)==3
            FWHM  = vertcat(FWHM1,currentLine(3));
        end
        
        currentLine = fgetl(fid);
    end
    
    fclose(fid);    
    % convert depth from cm to mm;
    ddd(i).depths = depth*10; 
    % save Z in MeV cm2 / g
    ddd(i).Z   = ionization;
    % extract peak position
    [~,idx]=max(ionization);
    ddd(i).peakPos = ddd(i).depths(idx);
    % add offset
    offset = -2.89;
    ddd(i).offset = ones(length(ionization),1)*offset;
    
    %% calculate sigmas and weights used to model lateral profile
    if ~isempty(FWHM1)
        % convert full width half maximum to sigma
        sigma1 = abs(FWHM1)/(2*sqrt(2*log(2)));
        sigma2 = abs(FWHM2)/(2*sqrt(2*log(2)));
        % according to Julian Streitz Thesis; 
        % add inital beam width to sigmas
        ddd(i).sigma1 = sqrt(sigmaSISsq(i,1) + sigma1(:,1).^2);                                                                
        ddd(i).sigma2 = sqrt(sigmaSISsq(i,1) + sigma2(:,1).^2);
        ddd(i).weight = weight; 
    elseif ~isempty(FWHM)
        sigma = abs(FWHM)/(2*sqrt(2*log(2)));
        ddd(i).sigma = sqrt(sigmaSISsq(i,1) + sigma(:,1).^2);     
    else           
        %% TODO: add analytical calculation of sigma according 
        %to Hong or the Highland formula
        ddd(i).sigma = ones(size(ddd(i).depths,1),1);
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
            Idx = 150;
        case {'p','H'}
            load protonBaseData   
            Idx = 250;
    end
    
    vIdx = round(linspace(10,length(baseData)-10,4));
    figure,
    for i = 1:length(vIdx)
        Energy = baseData(vIdx(i)).energy;
        [~,IdxDDD]=min(abs([ddd.energy]-Energy));
        subplot(2,2,i),plot(ddd(IdxDDD).depths,ddd(IdxDDD).Z,'r'),hold on
                       plot(baseData(vIdx(i)).depths,baseData(vIdx(i)).Z,'b');
                       legend({'new ddd - HIT','existing baseData'});
                       xlabel(' depth in mm '), ylabel('Z in MeVcm^2 / g'), title(['Energy: ' num2str(Energy)])
    end
    
    % plot sigmas againts depth
    figure, hold on 
    plot(ddd(Idx).depths,ddd(Idx).sigma1,'LineWidth',3)
    plot(ddd(Idx).depths,ddd(Idx).sigma2,'LineWidth',3)
    sigmaMix = (1-ddd(Idx).weight).*ddd(Idx).sigma1 + ddd(Idx).weight.*ddd(Idx).sigma2;
    plot(ddd(Idx).depths,sigmaMix,'LineWidth',3)
    grid on, grid minor, xlabel('depth in mm'),ylabel('sigma'),
    legend({'sigma1','sigma2','sigmas weighted'})
    title(['sigmas of a ' Identifier ' beam with energy ' num2str(ddd(Idx).energy) ' MeVcm^2/g'])
    
end