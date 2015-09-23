function [ deltaElectronRange ] = LEM_maxElectronRange(Energy,visBool)

% Energy in Mev/u

Identifier = 'OldScholz';

switch Identifier
    case 'OldScholz'
        deltaElectronRange = 0.05 * (Energy^1.7);
        
    case 'newScholz'
        deltaElectronRange = 0.062 * (Energy^1.7);
        
    case 'Geiss'
        deltaElectronRange = 0.04 * (Energy^1.5);
end



if visBool
   
    vEnergy = logspace(-1,3);
    vDeltaElectronRange = 0.05 * vEnergy.^1.7;
    figure,plot(vEnergy,vDeltaElectronRange),grid on, grid minor,
    set(gca,'Xscale','log'),set(gca,'Yscale','log'),
    xlabel('Energy in MeV/u'),ylabel('range of \delta electrons in µm = track diameter')
    title('range energy relationship of \delta electrons'), set(gca,'FontSize',14);
end


end

