function [ machine  ] = matRad_readBeamWidthHIT(machine,PathToXMLFile)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpLateralInfo script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ machine ] = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,FocusIdx,visBool)
%
% input
%   machine:           base data
%   pathTRiP:          path to TRiP folder for parsing the inital beam width
%   pathToSparseData:  path to sparse lateral double gauss data
%   Identifier:        either 'p','C','O' for parsing the correct beam
%                      width
%   focusIdx:          focus index (1-6) determines the initial beam width which will be added to the lateral
%                      sigma(s). If focusIdx is set to 0 no initial beam width will be added
%   visBool:           boolean if plots should be displayed or not - if visBool
%                      is 1, then press a key to continue with the next energy
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
% this call takes about 30
data = xml2struct(PathToXMLFile);

Fac = 2*sqrt(2*log(2));

% only this specific field is needed - meta information is omitted
GaussData = data.SingleGaussianBeamWidthSet.SingleGaussianBeamWidthData.SingleGaussianBeamWidth;

EnergyCnt = 0;
focusIx   = 0;
EnergyRef = 0;

for i = 1:length(GaussData)
   
        Key   = GaussData{1,i}.Key;
        Value = GaussData{1,i}.Value;
        
        Energy    = str2double(Key.Key1.Text);
        
        if Energy ~= EnergyRef
            focusIx = 1;
            EnergyCnt = EnergyCnt + 1;
            EnergyRef = Energy;
        else
            focusIx = focusIx + 1;
        end
        
        % just for cross checking
        SigmaInI  = str2double(Key.Key2.Text);
        FWHMini   = SigmaInI*Fac;
        SigmaInI2 = str2double(Key.Key3.Text);
        
        Cnt = 1;
        
        if abs(machine.data(EnergyCnt).energy - Energy) > 0.02
            warning('check if data is correctly added: discrepency in energies');
        end
        
        for k = 1:length(Value.BeamWidthX.Entry)           
   
                machine.data(EnergyCnt).initFocus(focusIx).dist(Cnt)  = str2double(Value.BeamWidthX.Entry{1,k}.x.Text);
                machine.data(EnergyCnt).initFocus(focusIx).sigma(Cnt) = str2double(Value.BeamWidthX.Entry{1,k}.y.Text);
              
                Cnt = Cnt + 1;
        end
        
        machine.data(EnergyCnt).initFocus(focusIx).unit   = 'mm';
        
   
end






