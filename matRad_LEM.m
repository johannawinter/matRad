clc,
clear
close all
addpath('C:\Users\wieserh\Documents\matRad\baseDataHIT');
load('C:\Users\wieserh\Documents\matRad\baseDataHIT\initialRBE.mat');
load('C:\Users\wieserh\Documents\matRad\baseDataHIT\dEdx.mat');
load('C12spc350.mat');
%% specify meta data
Nevent = 100; %1000, 10 000
Nruns = 500;
CellType = 1;
DepthStep = 50;

Smax = initialRBE(CellType).alpha+(2*initialRBE(CellType).beta*initialRBE(CellType).cut);
sParticle={'H','He','Li','Be','B','C'};
sParticleLong={'hydrogen','helium','lithium','beryllium','bor','carbon'};
SPCdat = C12spc350(DepthStep);
fluence = 0;
fluence_Z = zeros(6,1);
%% the total fluence = sum of the fluences of each particle species
for i = 1:length(sParticle)
    fluence_Z(i) = sum(SPCdat.(sParticle{i}).N);
    fluence = fluence +fluence_Z(i);
    
    %fluence = sum(SPCdat.(sParticle{i}).dNdE .* SPCdat.(sParticle{i}).dE);
end

%% calculate the number of hits on a cell nucleus
Anuc = pi*(initialRBE(CellType).rnucleus^2); %µm^2
Anuc = Anuc/(10000^2);  %cm^2
N_avg_hit = fluence*Anuc;

N_hit = poissrnd(2,1,Nevent);
%% convert the fluence into a probability distribution
PropDist = fluence_Z/fluence;

%% create a set of particles (particle type T(k), energy E(k))
Dabs = zeros(Nevent,1);
Nlethal = zeros(Nevent,1);
idx = 1:length(PropDist);

for j = 1:Nruns
   for i = 1:Nevent
       % first case: no particle hits a cell nucleus
       if N_hit(i) == 0
           
           Dabs(i) = 0;
           Nlethal(i) = 0;
           
       % second case: at least one particle hits a cell nucleus
       else
           
           %sampled n = N.Hit times to obtain a set (T(k), E(k))
           IdxParticle = datasample(idx,N_hit(i),'Replace',true,'Weights',PropDist);
           vDose = zeros(N_hit(i)+1,1);
           vDose(1)=0;
           vDoseTot = 0;
           % compute effective damage
           vNlet = zeros(N_hit(i)+1,1);
           vNletTot = 0;
           
           %sample the corresponding energy from the spectral distribution
           for m = 1:length(IdxParticle)
               
               energy=datasample(SPCdat.(sParticle{IdxParticle(m)}).Emid,1,'Replace',true);
                % compute absorbed dose and effective damage
                C = 1.602e-10/Anuc;
                dEdx_sampled = interp1(dEdx.(sParticle{IdxParticle(m)}).energy,dEdx.(sParticle{IdxParticle(m)}).dEdx,energy);
                vDose(m+1) = (dEdx_sampled*C);
              
                vDoseTot = vDoseTot + vDose(m+1);

                alphaIon = interp1([initialRBE(CellType).(sParticleLong{IdxParticle(m)}){1,2}.Energy],...
                                       [initialRBE(CellType).(sParticleLong{IdxParticle(m)}){1,2}.RBE],...
                                       dEdx_sampled,'linear','extrap')*initialRBE(CellType).alpha;

                Sk = alphaIon + (Smax-alphaIon)*(vDose(m)/initialRBE(CellType).cut);   

                   if m > 1
                       if vDose(m-1)<initialRBE(CellType).cut
                          Sk=Smax; 
                       end
                   end
               
                 vNlet(m+1) = Sk*dEdx_sampled*C;
                 vNletTot = vNletTot + vNlet(m+1);
                          
           end
           
           
       end
       
          Dabs(i) = vDoseTot;
          Nlethal(i) = vNletTot;
       
   end
   DabsAvg = sum(Dabs)/Nevent; 
   NlethalAvg = sum(Nlethal)/Nevent;
   
   
end



