clc, clear, close all

load('protonBasedataHIT');

eV2J = 1.60217662e-19;
f    =  1.60217662e-10;
offset = -2.89; % in mm

for i = 1:length(baseData)   
    
    xi = offset;
    Zi = interp1(baseData(i).depths,baseData(i).Z,xi);
    %Z
    [NewDepths,I] = sort([baseData(i).depths; xi]);
    NewZ = [baseData(i).Z; Zi];
    baseData(i).Z = NewZ(I)';
    
     %sigma1
     Zi = interp1(baseData(i).depths,baseData(i).sigma1,xi);
     NewZ = [baseData(i).sigma1; Zi];
     baseData(i).sigma1 = NewZ(I);
     
     %sigma2
     Zi = interp1(baseData(i).depths,baseData(i).sigma2,xi);
     NewZ = [baseData(i).sigma2; Zi];
     baseData(i).sigma2 = NewZ(I);
     
     %weight
     Zi = interp1(baseData(i).depths,baseData(i).weight,xi);
     NewZ = [baseData(i).weight; Zi];
     baseData(i).weight = NewZ(I);
     
     %add depth
     baseData(i).depths = NewDepths; 
    
     baseData(i).MeVmm2g= baseData(i).Z.*100;
     baseData(i).GyperParticle = baseData(i).MeVmm2g*f;
     baseData(i).GyperMilParticle = baseData(i).GyperParticle*1e6;
     baseData(i).Z = baseData(i).GyperMilParticle;
     baseData(i).offset = offset;    
end

Idx = 80;
figure,
subplot(131),plot(baseData(Idx).depths,baseData(Idx).MeVmm2g,'LineWidth',3),xlabel('depth'),ylabel('MeV mm^2/g')
grid on, grid minor,set(gca,'FontSize',14),title(['ddd of ' num2str(baseData(Idx).energy) 'MeV beam'])
subplot(132),plot(baseData(Idx).depths,baseData(Idx).GyperParticle,'LineWidth',3),xlabel('depth'),ylabel('Gy per primary')
grid on, grid minor,set(gca,'FontSize',14),title(['ddd of ' num2str(baseData(Idx).energy) 'MeV beam'])
subplot(133),plot(baseData(Idx).depths,baseData(Idx).Z,'LineWidth',3),xlabel('depth'),ylabel('Gy per million particles')
grid on, grid minor,set(gca,'FontSize',14),title(['ddd of ' num2str(baseData(Idx).energy) 'MeV beam'])


baseData = rmfield(baseData,'MeVmm2g');
baseData = rmfield(baseData,'GyperParticle');
baseData = rmfield(baseData,'GyperMilParticle');



%%
energyIx = 100;
depths = baseData(energyIx).depths + baseData(energyIx).offset;
Idx = depths >= -3;


figure,

plot(baseData(energyIx).depths,baseData(energyIx).Z,'Linewidth',2), hold on
plot(baseData(energyIx).depths,baseData(energyIx).Z,'bx','Linewidth',3)

plot(depths(Idx),baseData(energyIx).Z(Idx),'r','Linewidth',2)
plot(depths(Idx),baseData(energyIx).Z(Idx),'rx','Linewidth',3)


grid on, grid minor




