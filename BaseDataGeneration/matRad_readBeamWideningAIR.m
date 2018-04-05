function [ machine  ] = matRad_readBeamWideningAIR(machine,PathToXMLFile)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_readBeamWideningAIR script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ machine  ] = matRad_readBeamWideningAIR(machine,PathToXMLFile)
%
% input
%   machine:           base data file
%   PathToXMLFile:     path to LPD.xml which describes the beam widening in
%                      air. These files have been provided for protons and
%                      carbon ions from Katja.P. An additional field called
%                      initFocus containing an look up table will be added
%                      to the machine.data struct
%
% output
%   machine:           machine containing interpolated lateral double
%                      gaussian information
%   
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse all tags of the xml file

% this call takes about 30-40 seconds
data = xml2struct(PathToXMLFile);

FWHMtoSIGMA = 2*sqrt(2*log(2));

% only this specific field is needed - meta information is omitted
GaussData = data.SingleGaussianBeamWidthSet.SingleGaussianBeamWidthData.SingleGaussianBeamWidth;

EnergyCnt = 0;
focusIx   = 0;
EnergyRef = 0;

for i = 1:length(GaussData)
   
        Key    = GaussData{1,i}.Key;
        Value  = GaussData{1,i}.Value;
        Energy = str2double(Key.Key1.Text);
        
        if Energy ~= EnergyRef
            focusIx = 1;
            EnergyCnt = EnergyCnt + 1;
            EnergyRef = Energy;
        else
            focusIx = focusIx + 1;
        end
        
        %cross checking
        SigmaInI  = str2double(Key.Key2.Text);
        FWHMini   = SigmaInI*FWHMtoSIGMA;
        SigmaInIRef = str2double(Key.Key3.Text);
        if abs((SigmaInI/SigmaInIRef)-1) > 0.01
            warning('discrepancies detected');
        end
        Cnt = 1;
        
        if abs(machine.data(EnergyCnt).energy - Energy) > 0.02
            warning('check if data is correctly added: discrepency in energies');
        end
        
       
        for k = 1:length(Value.BeamWidthX.Entry)           
   
                machine.data(EnergyCnt).initFocus.dist(focusIx,Cnt)  = machine.meta.SAD - machine.meta.BAMStoIsoDist + str2double(Value.BeamWidthX.Entry{1,k}.x.Text);
                machine.data(EnergyCnt).initFocus.sigma(focusIx,Cnt) = str2double(Value.BeamWidthX.Entry{1,k}.y.Text);
              
                Cnt = Cnt + 1;
        end
        
        machine.data(EnergyCnt).initFocus.unit   = 'mm';
        
end






