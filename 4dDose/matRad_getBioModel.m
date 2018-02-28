function [ pln ] = matRad_getBioModel(pln)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% check for supported optimization methods using PHOTONS

if isequal(pln.radiationMode,'photons')
    
    if isequal(pln.propOpt.bioOptimization,'physical')
        bioParam.type             = 'physical';
    else
         warning(['matRad: Invalid biological optimization: ' pln.propOpt.bioOptimization ' for ' pln.radiationMode '; using physical optimization instead'])
    end

bioParam.radiationMode = pln.radiationMode;

%% check for supported optimization methods using PROTONS    
elseif isequal(pln.radiationMode,'protons')
    
    if isequal(pln.propOpt.bioOptimization,'none')
        bioParam.bioOpt             = false;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'none';
        bioParam.quantity           = 'physicalDose';
        
    elseif isequal(pln.propOpt.bioOptimization,'const_RBExD') 
        bioParam.bioOpt             = false;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'constant RBE as used clinically';
        bioParam.quantity           = 'RBExD';
        bioParam.constRBE           = 1.1;
        
    elseif isequal(pln.propOpt.bioOptimization,'LSM_effect')  
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'LSM';
        bioParam.description        = 'linear scaling model whereas the proton alpha is a linear function of the LETd; proton beta is constant';
        bioParam.quantity           = 'effect';
        bioParam.lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (fitted only for head and neck patients)
        bioParam.corrFacEntranceRBE = 0.5;   %[kev/mum]
        bioParam.upperLETThreshold  = 30;    %[kev/mum]
        bioParam.lowerLETThreshold  = 0.3;   %[kev/mum]
        
     elseif isequal(pln.propOpt.bioOptimization,'LSM_RBExD') 
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'LSM';
        bioParam.description        = 'linear scaling model whereas the proton alpha is a linear function of the LETd; proton beta is constant';
        bioParam.quantity           = 'RBExD';
        bioParam.lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (fitted for head and neck patients)
        bioParam.corrFacEntranceRBE = 0.5;   %[kev/mum]
        bioParam.upperLETThreshold  = 30;    %[kev/mum]
        bioParam.lowerLETThreshold  = 0.3;   %[kev/mum]
        
      elseif isequal(pln.propOpt.bioOptimization,'MCN_effect') 
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'MCN';
        bioParam.description        = 'a phenomenological relative biological effectiveness (RBE) model for proton therapy based on all published in vitro cell survival data';
        bioParam.quantity           = 'effect';
        bioParam.p0                 = 0.999064; % according to https://www.ncbi.nlm.nih.gov/pubmed/26459756
        bioParam.p1                 = 0.35605;
        bioParam.p2                 = 1.1012;
        bioParam.p3                 = -0.0038703;  
        
      elseif isequal(pln.propOpt.bioOptimization,'MCN_RBExD') 
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'MCN';
        bioParam.description        = 'a phenomenological relative biological effectiveness (RBE) model for proton therapy based on all published in vitro cell survival data';
        bioParam.quantity           = 'RBExD';
        bioParam.p0                 = 0.999064; % according to https://www.ncbi.nlm.nih.gov/pubmed/26459756
        bioParam.p1                 = 0.35605;
        bioParam.p2                 = 1.1012;
        bioParam.p3                 = -0.0038703;
        
       elseif isequal(pln.propOpt.bioOptimization,'WED_effect') 
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'WED';
        bioParam.description        = 'A model for the relative biological effectiveness of protons: the tissue specific parameter alpha/beta of photons is a predictor for the sensitivity to LET changes.';
        bioParam.quantity           = 'effect';
        bioParam.p0                 = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
        bioParam.p1                 = 0.434;
        bioParam.p2                 = 1;
        bioParam.p3                 = 0;
        
       elseif isequal(pln.propOpt.bioOptimization,'WED_effect') 
        bioParam.bioOpt             = true;
        bioParam.type               = pln.propOpt.bioOptimization;
        bioParam.model              = 'WED';
        bioParam.description        = 'A model for the relative biological effectiveness of protons: the tissue specific parameter alpha/beta of photons is a predictor for the sensitivity to LET changes.';
        bioParam.quantity           = 'RBExD';
        bioParam.p0                 = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
        bioParam.p1                 = 0.434;
        bioParam.p2                 = 1;
        bioParam.p3                 = 0; 
        
    else
         warning(['matRad: Invalid biological optimization: ' pln.propOpt.bioOptimization ' for ' pln.propOpt.radiationMode '; using const_RBExD optimization instead'])
         % back up solution
        bioParam.bioOpt           = true;
        bioParam.type             = 'const_RBExD';
        bioParam.model            = 'constantRBE';
        bioParam.quantity         = 'RBExD';
        bioParam.constRBE         = 1.1;
    end
    
%% check for supported optimization methods using CARBON ions  
elseif isequal(pln.radiationMode,'carbon')
    
    if isequal(pln.propOpt.bioOptimization,'physical')
        bioParam.bioOpt           = false;
        bioParam.type             = 'physical';
        bioParam.model            = 'none';
        bioParam.quantity         = 'physicalDose';
        
    elseif isequal(pln.propOpt.bioOptimization,'LEMIV_effect')  
        bioParam.bioOpt           = true;
        bioParam.type             = pln.propOpt.bioOptimization;
        bioParam.model            = 'LEMIV';
        bioParam.quantity         = 'effect';
        bioParam.BeamMixingModel  = 'ZaiderRossi';
        bioParam.dummy = 1;
        
    elseif isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD') 
        bioParam.bioOpt           = true;
        bioParam.type             = pln.propOpt.bioOptimization;
        bioParam.model            = 'LEMIV';
        bioParam.quantity         = 'RBExD';
        bioParam.BeamMixingModel  = 'ZaiderRossi';
        
    else
         warning(['matRad: Invalid biological optimization: ' pln.propOpt.bioOptimization ' for ' pln.radiationMode '; using LEMIV_RBExD optimization instead'])
        bioParam.bioOpt           = true;
        bioParam.type             = pln.propOpt.bioOptimization;
        bioParam.Model            = 'LEMIV';
        bioParam.quantity         = 'RBExD';
        bioParam.BeamMixingModel  = 'ZaiderRossi';
    end


end


pln.bioParam = bioParam; 
