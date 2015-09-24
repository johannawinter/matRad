function baseData = matRad_getDDDfromTxt(Identifier,basePath,FocusIdx,Offset,visBool)
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
    [Sigma_SISsq,vEnergySIS] = matRad_getSigmaSISsq(Identifier,basePath,FocusIdx);
catch
    warning('Could not read sis files containing the inital beam widths');
    Sigma_SISsq = 0;
end

% check the content of the directory
Files = dir([basePath relPath filesep '*.ddd']);

% parse all ddd files
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
    baseData(i).energy = strread(currentline,'!energy %f');

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
    
    % read the whole file
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
    baseData(i).depths = depth*10; 
    % save Z in MeV cm2 / g
    baseData(i).Z   = ionization;
    % extract peak position
    [~,idx]=max(ionization);
    baseData(i).peakPos = baseData(i).depths(idx);
    
    baseData(i).FWHM1 =FWHM1;
    baseData(i).FWHM2 =FWHM2;
    baseData(i).weight =weight;
   
end   

% sort content according to energy - in ascending order;
[~,IdxEnergy] = sort([baseData.energy]);
baseData=baseData(IdxEnergy);
%add offset to each entry
for i = 1:length(baseData)
   baseData(i).offset = Offset; 
end

if sum(abs([baseData.energy]'-vEnergySIS))>1
   disp('check if foki and ddd exist on the same energy levels');
end

% calculate sigmas and weights used to model lateral profile
for i = 1:length(baseData)
    
    if ~isempty(FWHM1)
        % convert full width half maximum to sigma
        sigma1 = baseData(i).FWHM1/(2*sqrt(2*log(2)));
        sigma2 = baseData(i).FWHM2/(2*sqrt(2*log(2)));
        % add inital beam width to sigmas
        baseData(i).sigma1 = sqrt(Sigma_SISsq(i,1) + sigma1(:,1).^2);                                                                
        baseData(i).sigma2 = sqrt(Sigma_SISsq(i,1) + sigma2(:,1).^2);
        
    elseif ~isempty(FWHM)
        
        sigma = abs(FWHM)/(2*sqrt(2*log(2)));
        baseData(i).sigma = sqrt(Sigma_SISsq(i,1) + sigma(:,1).^2);   
        
    else           
        %% TODO: add analytical calculation of sigma according 
        %to Hong or the Highland formula
        baseData(i).sigma = ones(size(baseData(i).depths,1),1);
    end
end

% remove FWHM fields
baseData = rmfield(baseData,'FWHM1');
baseData = rmfield(baseData,'FWHM2');




%% plot random baseData 
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
        [~,IdxDDD]=min(abs([baseData.energy]-Energy));
        subplot(2,2,i),plot(baseData(IdxDDD).depths,baseData(IdxDDD).Z,'r'),hold on
                       plot(baseData(vIdx(i)).depths,baseData(vIdx(i)).Z,'b');
                       legend({'new ddd - HIT','existing baseData'});
                       xlabel(' depth in mm '), ylabel('Z in MeVcm^2 / g'), title(['Energy: ' num2str(Energy)])
    end
    
    % plot sigmas againts depth
    figure, hold on 
    plot(baseData(Idx).depths,baseData(Idx).sigma1,'LineWidth',3)
    plot(baseData(Idx).depths,baseData(Idx).sigma2,'LineWidth',3)
    sigmaMix = (1-baseData(Idx).weight).*baseData(Idx).sigma1 + baseData(Idx).weight.*baseData(Idx).sigma2;
    plot(baseData(Idx).depths,sigmaMix,'LineWidth',3)
    grid on, grid minor, xlabel('depth in mm'),ylabel('sigma'),
    legend({'sigma1','sigma2','sigmas weighted'})
    title(['sigmas of a ' Identifier ' beam with energy ' num2str(baseData(Idx).energy) ' MeVcm^2/g'])
    
end