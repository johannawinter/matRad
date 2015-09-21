clc, clear, close all
addpath('C:\Users\wieserh\Documents\matRad');
load('protonBaseDataHIT');

eVtoJ   = 1.60217662e-19;
Factor  = eVtoJ * 1e9;


offset = -2.89; % in mm

for i = 1:length(baseData)   
%      baseData(i).MeVmm2g= baseData(i).Z.*100;
%      baseData(i).GyperParticle = baseData(i).MeVmm2g*Factor;
%      baseData(i).GyperMilParticle = baseData(i).GyperParticle*1e6;
%      baseData(i).Z = baseData(i).GyperMilParticle;
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




