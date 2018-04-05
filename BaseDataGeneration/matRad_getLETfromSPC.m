function [ machine ] = matRad_getLETfromSPC( machine,Identifier,pathTRiP )


BoolSaveSPCasMAT = true;

switch Identifier
    case {'c','C'}
        relPath   = [filesep '12C' filesep 'RF3MM' filesep];
        Particles = {'H','He','Li','Be','B','C'};
    case {'p','P','h''H'}
        relPath = [filesep 'p' filesep' 'RF0MM' filesep ];
        Particles = {'H'};
    case {'o','O'} 
         error('unkown particle type')
    otherwise
        error('unkown particle type')
end

pathToSPCFiles = [pathTRiP filesep 'SPC' relPath];

%% convert spc to mat files
DirInfo = dir([pathToSPCFiles '*.spc']);
for i = 1:length(DirInfo)
    
 [metaSPC, SPCdata] = matRad_readSPC([pathToSPCFiles filesep DirInfo(i).name]);
 sSPC(i).energy = metaSPC.energy;
 sSPC(i).data   = SPCdata; 
 sSPC(i).meta   = metaSPC;

 if BoolSaveSPCasMAT
     SPC.meta = metaSPC;
     SPC.data = SPCdata;
     fileName = DirInfo(i).name;
     save([pathToSPCFiles filesep fileName(1:end-4) '.mat'],'SPC');
 end
 
end

%% load dEdx and calculate LET
[MetadEdx, dEdx ] = matRad_readdEdx(pathTRiP);

LET_MODE = 'dose_averaged';

for ixEnergy = 1:length(sSPC)
    
    vDepth      = [sSPC(ixEnergy).data.depths];
    vLET_Accum  = zeros(length(vDepth),1);
    vDose_Accum = zeros(length(vDepth),1); 
    
    for ixPart = 1:length(Particles)
        
        vLET_Numerator  = zeros(length(vDepth),1); 
        vDose          = zeros(length(vDepth),1); 
        
         for ixDepth = 1:size(sSPC(ixEnergy).data,2)
            
              ix = find(sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).Emid >= dEdx.(Particles{ixPart}).Energy(1) & ...
                        sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).Emid < dEdx.(Particles{ixPart}).Energy(end));
              
              dEdxInterp = interp1(dEdx.(Particles{ixPart}).Energy,...
                                   dEdx.(Particles{ixPart}).dEdx,...
                                   sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).Emid(ix),'linear')'; 
              switch   LET_MODE
                  case {'dose_averaged'}
                          vLET_Numerator(ixDepth) = (sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).N(ix) *(dEdxInterp).^2);
                          vDose(ixDepth)          = (sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).N(ix) *(dEdxInterp));
                  case {'track_averaged'}
                         
                      vLET_Numerator(ixDepth) = (sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).N(ix) *(dEdxInterp));
                          vDose(ixDepth)          = sum(sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).N(ix));
                          
%                            vLET_Numerator(ixDepth) = (sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).dNdE(ix).*...
%                                                          sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).dE(ix)) *(dEdxInterp).^2;
%                            vDose(ixDepth)          = (sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).dNdE(ix).*...
%                                                          sSPC(ixEnergy).data(ixDepth).(Particles{ixPart}).dE(ix)) *dEdxInterp;
                  otherwise
                      error('mode doesnt exitst')
              end
         end
         
         vLET_Accum  =  vLET_Accum  + vLET_Numerator;
         vDose_Accum =  vDose_Accum + vDose;
    end
    
    sSPC(ixEnergy).LET = (vLET_Accum./(vDose_Accum*10))'; % convert from MeV/cm to kEV/µm
end


%% interpoalte LET 
% figure,
% for ixEnergy = 1:length(sSPC)
%     plot([sSPC(ixEnergy).data.depths],sSPC(ixEnergy).LET),hold on
% end
% set(gca,'YScale','log')

for idxE = 1:length(machine.data)
    
   [~,idxSorted] = sort(abs(machine.data(idxE).energy - [sSPC.energy]));
   
   Idx = sort(idxSorted(1:3));
   LengthDepth = size(machine.data(idxE).depths,1);
   vLET = zeros(LengthDepth,3); 
   
   vLET(:,1) = interp1([sSPC(Idx(1)).data.depths]./sSPC(Idx(1)).data(1).peakPos,...
                        sSPC(Idx(1)).LET,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');
   
   vLET(:,2) = interp1([sSPC(Idx(2)).data.depths]./sSPC(Idx(2)).data(1).peakPos,...
                         sSPC(Idx(2)).LET,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vLET(:,3) = interp1([sSPC(Idx(3)).data.depths]./sSPC(Idx(3)).data(1).peakPos,...
                         sSPC(Idx(3)).LET,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   for k = 1:LengthDepth
       machine.data(idxE).LET(k,1) = interp1([sSPC(Idx).energy],vLET(k,:),machine.data(idxE).energy,'linear','extrap');
   end
     
end



end

