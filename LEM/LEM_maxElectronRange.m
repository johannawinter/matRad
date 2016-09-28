function [ deltaElectronRange ] = LEM_maxElectronRange(Energy,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_maxElectronRange
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the maximum delta elctron range

Identifier = 'OldScholz';

switch Identifier
    case 'OldScholz'
        ElectronRangeModel = @(vEnergy)(0.05 * (Energy^1.7));
        
    case 'newScholz'
        ElectronRangeModel = @(vEnergy)(0.062 * (Energy^1.7));
        
    case 'Geiss'
        ElectronRangeModel = @(vEnergy)(0.04 * (Energy^1.5));
end

deltaElectronRange = ElectronRangeModel(Energy);

if visBool 
    vEnergy = logspace(-1,3);
    vDeltaElectronRange = ElectronRangeModel(vEnergy);
    figure,plot(vEnergy,vDeltaElectronRange),grid on, grid minor,
    set(gca,'Xscale','log'),set(gca,'Yscale','log'),
    xlabel('Energy in MeV/u','Interpreter','Latex'),...
    ylabel('range of $\delta$ electrons in $\mu$m = track diameter','Interpreter','Latex')
    title('range energy relationship of $\delta$ electrons','Interpreter','Latex'), set(gca,'FontSize',14);
end


end

