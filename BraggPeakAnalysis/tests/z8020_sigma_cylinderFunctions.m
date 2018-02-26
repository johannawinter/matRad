% find connection between z8020 and sigma with parabolic cylinder function

clear
load('protons_HIT_APM')
p = 1.77;
alphaWater = 2.2e-3;
alphaPrimeWater = 0.087;

ixBeamEnergy = 70;

R0 = alphaWater * machine.data(ixBeamEnergy).energy ^p;   % [cm]
z = machine.data(ixBeamEnergy).depths;                    % [mm]
sigmaSq = alphaPrimeWater * (p^3*alphaWater^(2/p))/(3*p-2) * R0.^(3-2/p);
sigma = sqrt(sigmaSq) *10;          % [mm]

chi = -(R0*10-z)/sigma;             % [mm]
relDose = exp(-chi(48:104)/4) .* cylFun(-1/p,-chi(48:104));

figure
hold on
plot(z(48:104),relDose)
