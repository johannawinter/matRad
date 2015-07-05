function [ sData ] = matRad_ParseBioDataHIT(PathToHITBaseData,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad fucntion to parse HIT dEdx*alpha and dEdx*sqrtBeta curves from
% spc data files
%
% call
%   sData = matRad_ParseBioDataHIT(PathToBaseData,visBool)
% example
%   sData = matRad_ParseBioDataHIT('C:\Users\admin\baseData',visBool)
%   e.g. the folder baseData contains two subfolders AB2 and AB10.
%
% input
%   PathToHITBaseData:    path to folder containing dEdx,initial RBE and
%                         spc files in *.mat format
%   visBool:              toggle on/off visualization (optional)
%
% output
%   sData:                 returns a cell array, whereat each cell
%                          represents the bio data for one specific cell line 
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Spectra = {'hydrogen','helium','lithium','beryllium','bor','carbon','nitrogen','oxygen','fluor','neon'};
sParticles = {'H','He','Li','Be','B','C'};
sParticleLong = {'hydrogen','helium','lithium','beryllium','bor','carbon'};
sColor={'red','green','blue','red','green','blue','black'};
sLineSpec={'--','--','--','-','-' ,'-' ,'-'};
sLineSpec2 = {':',':',':','-.','-.','-.'};
DimEnergy = 37;
CntEnergy = 1;

%% load data

load([PathToHITBaseData filesep 'dEdx.mat']);
load([PathToHITBaseData filesep 'initialRBE.mat']);
load([PathToHITBaseData filesep 'carbonBaseDataHIT.mat']);

%% get alphas and betas
path = [PathToHITBaseData filesep];
% loop over all cell lines
for j = 1:length(initialRBE)
 % loop over all energies
 for i=80:10:440
     
  if i/10<=9
        filename = ['C12spc' '0' num2str(i)];
  else
        filename = ['C12spc' num2str(i)];
  end
     
  load([path filename]);
  fName=fieldnames(SPC);
  SPC = SPC.(fName{1});
  vDepth = [SPC(:).depth];
  % extract meta data for current cell line
  RBE_ini = initialRBE(j);
  alpha_x = initialRBE(j).alpha;
  beta_x = initialRBE(j).beta;
  Dcut = initialRBE(j).cut;
  Smax = alpha_x+(2*beta_x)*Dcut;
  Anuc = pi*(initialRBE(j).rnucleus^2); %µm^2
  Anuc = Anuc/(10000^2);
  % allocate variables
  alpha_numeratorRapid = zeros(length(vDepth),1);
  beta_numeratorRapid = zeros(length(vDepth),1);
  dose_accum = zeros(length(vDepth),1);
  
  for k = 1:length(sParticles)
      
      RBE_ini_z = RBE_ini.(sParticleLong{k}){1,2};
      alpha_ion = ([RBE_ini_z.RBE].*alpha_x)';
      beta_ion = (Smax-alpha_ion)./(2*Dcut);
      % rapid calculation according to Krämer and Scholz 2006
      dEdx_interp_RBE = interp1(dEdx.(sParticles{k}).energy,dEdx.(sParticles{k}).dEdx,[initialRBE(j).(sParticleLong{k}){1,2}.Energy],'pchip','extrap');
      d1 = ((1.602189e-10 .* dEdx_interp_RBE )/ Anuc);
      % S1 is the surviving fraction for a single particle traversal
      S1 = exp(-alpha_ion'.*d1);
      alpha_ion_rapid = (1-S1)./d1;
      % calculate scaling factor
      f = alpha_ion_rapid./alpha_ion';
      beta_ion_rapid=(f.^2).*beta_ion';
      %initialize some vectors
      alpha_Z_rapid = zeros(length(vDepth),1);
      beta_Z_rapid = zeros(length(vDepth),1);
      dose_Z = zeros(length(vDepth),1);
        
      
      for depth = 1:length(vDepth); 
        
        Fluence = (SPC(depth,1).(sParticles{k}).N);
        SP_interp = interp1(dEdx.(sParticles{k}).energy,dEdx.(sParticles{k}).dEdx,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap')';
        dose_Z(depth)= (Fluence*SP_interp);            
        
        alpha_ion_rapid_interp=interp1([RBE_ini_z.Energy],alpha_ion_rapid,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap');
        beta_ion_rapid_interp=interp1([RBE_ini_z.Energy],beta_ion_rapid,SPC(depth,1).(sParticles{k}).Emid,'pchip','extrap');
        
        alpha_Z_rapid(depth)= (alpha_ion_rapid_interp.*SP_interp')*Fluence';
        beta_Z_rapid(depth)= real((sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence');
          
      
      end
      
      alpha_numeratorRapid = alpha_numeratorRapid +alpha_Z_rapid;
      beta_numeratorRapid  = beta_numeratorRapid+ beta_Z_rapid;
      dose_accum = dose_accum + dose_Z;
      
  end
     
  sData{1,j}{CntEnergy}.energy = i;
  sData{1,j}{CntEnergy}.alphaBetaRatio = alpha_x/beta_x;
  sData{1,j}{CntEnergy}.depths = (vDepth')*10;
  sData{1,j}{CntEnergy}.alpha = alpha_numeratorRapid./dose_accum;
  sData{1,j}{CntEnergy}.beta = (beta_numeratorRapid./dose_accum).^2;
  sData{1,j}{CntEnergy}.alphaX = alpha_x;
  sData{1,j}{CntEnergy}.betaX = beta_x;
 
  %% check if division causes peaks in the tail
%   figure, hold on
%   plot(vDepth,sData{1,j}{CntEnergy}.beta)
%   [maxVal,peakPos] = max(dose_accum);
%   SecurityBuffer = 4;
%   Minimum = sData{1,j}{CntEnergy}.alpha(peakPos+SecurityBuffer);
%   for CntDepth = (peakPos+SecurityBuffer):1:length(vDepth)
%       if Minimum < sData{1,j}{CntEnergy}.alpha(CntDepth)
%         sData{1,j}{CntEnergy}.alpha(CntDepth) = Minimum;
%       else
%           Minimum = sData{1,j}{CntEnergy}.alpha(CntDepth);
%       end
%   end
   % plot(vDepth,sData{1,j}{CntEnergy}.alpha,'r')
  
    
    CntEnergy = CntEnergy+1;
 end
 CntEnergy = 1;
 
 
end

% convert cell to struct
for i = 1:length(sData)
 sData{1,i}= (cell2mat(sData{1,i}));
end


end
