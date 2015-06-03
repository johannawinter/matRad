clc,
clear all
close all
addpath('C:\Users\wieserh\Documents\matRad\baseDataHIT');
load('initialRBE');
load('dEdx');
load('C12spc350.mat');
%% specify meta data
Nevent = 100; %1000, 10 000
Nruns = 500;
CellType = 1;
DepthStep = 50;
sParticle={'H','He','Li','Be','B','C'};
SPCdat = C12spc350(DepthStep);
fluence = 0;
%% the total fluence = sum of the fluences of each particle species
for i = 1:length(sParticle)
    fluence = sum(SPCdat.(sParticle{i}).dNdE .* SPCdat.(sParticle{i}).dE);
    fluence2 = sum(SPCdat.(sParticle{i}).N);
end
fluence = fluence * 100;
%% calculate the number of hits on a cell nucleus
Anuc = pi*(initialRBE(CellType).rnucleus^2); %µm^2
Anuc = Anuc/(10000^2);  %cm^2
N_avg_hit = fluence*Anuc;

N_hit = poissrnd(2,1,Nevent);
%% convert the fluence into a probability distribution
PropDist = SPCdat.(sParticle{i}).N/sum(SPCdat.(sParticle{i}).N);

%% create a set of particles (particle type T(k), energy E(k))
Dabs = zeros(Nevent,1);
Nlethal = zeros(Nevent,1);

idx = 1:lenght(ProbDist);

for j = 1:Nruns
   for i = 1:Nevent
       % first case: no particle hits a cell nucleus
       if N_hit(i) == 0
           Dabs(i) = 0;
           Nlethal(i) = 0;
           
       % second case: at least one particle hits a cell nucleus
       else
           
           %sampled n = N.Hit times to obtain a set (T(k), E(k))
           datasample(data,N_hit(i),
           
       end
       
       
       
   end
    
end



