%% Comparison of shape of MC simulation data: no phantom vs. phantom A

load('RicData.mat');
offset_1 = 0.833;
offset_2 = 1.1;
offset_3 = 3.05;
Lz_A = 30;

coords_pE700sim = 10*mean(pE700sim(:,[1 2]),2);
coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

figure(10);
hold on; 
plot(coords_pE700sim - offset_2,pE700sim(:,3)./max(pE700sim(:,3)),'g+'); 
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3)./max(pE70Asim(:,3)),'b+');
axis([0 150 0 1]);
legend({'MC without phantom','MC phantom A'},'location','northwest');

% Shift phantom A curve so that the peaks are at the same position
[~,maxAIx] = max(pE70Asim(:,3));
[~,max0Ix] = max(pE700sim(:,3));
offsetPeak = coords_pE700sim(max0Ix) - coords_pE70Asim(maxAIx); % -20.4000
offsetTemp = offsetPeak + Lz_A;                                 % -21.5000

figure (11);
hold on;
plot(coords_pE700sim - offset_2, pE700sim(:,3)./max(pE700sim(:,3)),'g+'); 
plot(coords_pE70Asim - Lz_A - offset_2 + offsetTemp, pE70Asim(:,3)./max(pE70Asim(:,3)),'b+');
axis([0 150 0 1]);
legend({'MC without phantom','MC phantom A'},'location','northwest');
box on;


%% Comparison of different offset_1 on pE700
load('RicData.mat');
offset_1 = 0.833;
offset_2 = 1.1;
coords_pE700sim = 10*mean(pE700sim(:,[1 2]),2);
coords_pE700exp = 10*mean(pE700exp(:,[1 2]),2);

figure(40); 
hold on;
plot(coords_pE700exp-.815-offset_2,pE700exp(:,3)./max(pE700exp(:,3)),'b'); 
plot(coords_pE700exp-.833-offset_2,pE700exp(:,3)./max(pE700exp(:,3)),'r');
plot(coords_pE700sim-offset_2,pE700sim(:,3)./max(pE700sim(:,3)),'g');
axis([40 150 0 1]);
legend({'pE700 with offset\_1 = .815','pE700 with offset\_1 = .833','pE700sim'},'location','northwest');
box on;


%% Find offset_3 for processRicDataProton
load('RicData.mat');

[~,maxSIx] = max(pE700expStephan(:,2));
[~,max0Ix] = max(pE700exp(:,3));
offsetStephan = pE700expStephan(maxSIx,1) - coords_pE700exp(max0Ix) + offset_1 + offset_2;

figure(5);
hold on;
plot(coords_pE700exp-offset_1-offset_2,pE700exp(:,3)./max(pE700exp(:,3)),'bx');
plot(pE700expStephan(:,1)-offsetStephan,pE700expStephan(:,2)./max(pE700expStephan(:,2)),'r');
axis([0 150 0 1]);
legend({'exp data','Stephans data'},'location','northwest');


%% Comparison of curve shape of phantom A
load('RicData.mat');

offset_1 = 0.833;
offset_2 = 1.1;
offset_3 = 3.05;
Lz_A = 30;

coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);
coords_pE70Aexp = 10*mean(pE70Aexp(:,[1 2]),2);
coords_matRad = .5:1:350;

matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

[~,maxAIx] = max(pE70Aexp(:,3));
[~,maxMatRadA] = max(matRad_idd_pE70A);
offsetA = coords_pE70Aexp(maxAIx) - coords_matRad(maxMatRadA) - offset_1 - offset_2;

figure(36);
hold on;
plot(coords_pE70Aexp - offset_1 - offset_2,pE70Aexp(:,3)./max(pE70Aexp(:,3)),'bx');
plot(coords_matRad + offsetA,matRad_idd_pE70A./max(matRad_idd_pE70A),'r');
axis([0 150 0 1]);
legend({'pE70Aexp','matRad with shift'},'location','northwest');
