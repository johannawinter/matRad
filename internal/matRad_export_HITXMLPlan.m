%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Description: Function to export a matRad-based treatment plan in the HIT-XML
%               format based on the R-package HITXML from Steffen Greilich
%               available at https://r-forge.r-project.org/projects/hitxml/
%
%  Author:      Lucas Norberto Burigo
%
%  Please, report any problems to l.burigo@dkfz.de
%
%  Note1: the parameter minNbParticlesSpot is used to exclude spots with
%         number of particles below the threshold which can be delivered at HIT
%  Note2: the assigment of parameters for TxTable and BAMS needs to be verified
%  Note3: TxRoom fixed to Room1Fixed90
%  Note4: parameters for Beam, Patient and TxInitiation assigned to dumb values
%  Note5: for execution with Octave, check the Xerces Java Parser path below
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function XMLplan = matRad_export_HITXMLPlan(planFilename,minNbParticlesSpot)

if (nargin < 2)
  minNbParticlesSpot=0;
end

disp('HITXML exporter: Exporting plan in the HITXML format')

% get all data from workspace
pln       = evalin('base','pln');
stf       = evalin('base','stf');
resultGUI = evalin('base','resultGUI');

if ~strcmp(pln.radiationMode,'protons') && ~strcmp(pln.radiationMode,'carbon')
  error('HITXML plan for this radiationMode not supported!');
end

if ~strncmp(pln.machine,'HIT',3) 
  error('HITXML exporter only supports machine HIT');
end

% Compute bixel index as it is done on matRad_calParticleDose
counter = 0;

for i = 1:pln.numOfBeams; % loop over all beams

    for j = 1:stf(i).numOfRays % loop over all rays

        if ~isempty(stf(i).ray(j).energy)

            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                
                % remember beam and  bixel number
                tmp.beamNum(counter)  = i;
                tmp.rayNum(counter)   = j;
                tmp.bixelNum(counter) = k;
            end
        end
    end
end

% load machine file, required to read machine.data.initFocus.SisFWHMAtIso
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end
availableEnergies = [machine.data.energy];

% OCTAVE dependencies
if exist('OCTAVE_VERSION','builtin');

  pkg load io;

%  % check Xerces Java Parser dependencies
%  [status cmdout] = system('if test -e /usr/share/java/xerces-j2.jar; then exit 0; echo exit 2; fi');
%  if status == 2
%   error('Xerces Java Parser dependencies not found');
%    return;
%  end
%
%  [status cmdout] = system('if test -e /usr/share/java/xml-commons-apis.jar; then exit 0; echo exit 2; fi');
%  if status == 2
%   error('Xerces Java Parser dependencies not found');
%    return;
%  end
%
%  javaaddpath ('/usr/share/java/xerces-j2.jar');
%  javaaddpath ('/usr/share/java/xml-commons-apis.jar');

  if ~isempty(getenv('XERCES_DIR'))
    javaaddpath([getenv('XERCES_DIR'),'/xercesImpl.jar']);
    javaaddpath([getenv('XERCES_DIR'),'/xml-apis.jar']);
  else
    javaaddpath('/usr/local/sw/modules/source/xerces-java/xerces-2_11_0/xercesImpl.jar');
    javaaddpath('/usr/local/sw/modules/source/xerces-java/xerces-2_11_0/xml-apis.jar');
  end
end

for beamNb = 1:pln.numOfBeams % loop over beams
  XMLplan(beamNb).beamNb = beamNb;
  XMLplan(beamNb).couchAngle = stf(beamNb).couchAngle;
  XMLplan(beamNb).gantryAngle = stf(beamNb).gantryAngle;
  if strcmp(pln.radiationMode,'protons')
    XMLplan(beamNb).ion = '1H';
    XMLplan(beamNb).ionZ=1;
    XMLplan(beamNb).ionA=1;

  elseif strcmp(pln.radiationMode,'carbon')
    XMLplan(beamNb).ion = '12C';
    XMLplan(beamNb).ionZ=6;
    XMLplan(beamNb).ionA=12;
  end

  % helper function for energy selection
  round2 = @(a,b)round(a*10^b)/10^b;

  % Find IES values
  iesArray = [];
  for rayNb=1:stf(beamNb).numOfRays
    iesArray = unique([iesArray stf(beamNb).ray(rayNb).energy]);
  end

  iesNb=0; % used to get rid off IES in the iesArray with no particles (weigth = 0)
  for iesArrayIx=1:length(iesArray) % find focus and "Voxel (x,y,nbParticles)" for each IES
    clear iesFocus;
    clear iesEnergy;
    newIES=true;
    iesEnergy = iesArray(iesArrayIx);
    iesFocus = 0;

    for rayNb=1:stf(beamNb).numOfRays % loop over rays

      % find index of used energy (round to keV for numerical reasons
      bixelNb = find(round2(stf(beamNb).ray(rayNb).energy,4) == round2(iesEnergy,4)==1);
      energyIx = find(round2(iesEnergy,4) == round2(availableEnergies,4)==1);

      if length(bixelNb)==1 % one IES found

        bixelIndex = find([tmp.beamNum==beamNb & tmp.rayNum==rayNb & tmp.bixelNum==bixelNb]==1);

        voxel_nbParticles = resultGUI.w(bixelIndex);
        voxel_nbParticles = round(1e6*voxel_nbParticles);

        % check whether there are (enough) particles for beam delivery
        if (voxel_nbParticles>minNbParticlesSpot)

          rayPos_bev = stf(beamNb).ray(rayNb).rayPos_bev;
          % matRad, x-y(beam)-z
          % HITXML, x-y-z(beam)
          %
          %        -X
          %        |
          %  -Y -- o -- +Y
          %        |
          %        +X
          %
          voxel_x = -rayPos_bev(3);
          voxel_y = rayPos_bev(1);
  
          rayIESenergy=stf(beamNb).ray(rayNb).energy(bixelNb);
          rayIESfocusIx=stf(beamNb).ray(rayNb).focusIx(bixelNb);

          if newIES % new IES
            voxNb = 1;  
            iesNb = iesNb + 1;
            iesFocusIx = rayIESfocusIx;
            iesFocus = machine.data(energyIx).initFocus.SisFWHMAtIso(rayIESfocusIx);

%            XMLplan(beamNb).IES(iesNb).number = iesNb;
            XMLplan(beamNb).IES(iesNb).energy = iesEnergy;
            XMLplan(beamNb).IES(iesNb).focus = iesFocus;

            newIES=false;

            submachines{iesNb}.energyIx = energyIx;
            submachines{iesNb}.energy = iesEnergy;
            submachines{iesNb}.focusIx = iesFocusIx;
            submachines{iesNb}.focus = iesFocus;

          else % check whether the focus from the new bixel is the same as the current focus
            iesFocusIxNewBixel = rayIESfocusIx;
            if ( iesFocusIx ~= iesFocusIxNewBixel)
              warndlg('ATTENTION: bixels in same IES with different foci!');
              warndlg('ATTENTION: only one focus is kept for consistence with the Syngo TPS at HIT!');
            end
          end % new IES

          XMLplan(beamNb).IES(iesNb).voxel(voxNb).x = voxel_x;
          XMLplan(beamNb).IES(iesNb).voxel(voxNb).y = voxel_y;
          XMLplan(beamNb).IES(iesNb).voxel(voxNb).particles = voxel_nbParticles;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,1) = voxel_x;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,2) = voxel_y;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,3) = voxel_nbParticles;

          voxNb = voxNb + 1;
          
        end % there are particles for given voxel
      elseif length(bixelNb)>1
        error('Unexpected number of IES in the same ray.');
        return;
      end % one IES found
    end % loop over rays
  end % find focus and "Voxel (x,y,nbParticles)" for each IES


  filename = sprintf('PBP_%02d_%s.xml',beamNb-1,planFilename);

if exist('OCTAVE_VERSION','builtin');
  % OCTAVE
  docNode = javaObject ('org.apache.xerces.dom.DocumentImpl');
  docNode.appendChild (docNode.createElement ('PTTxPlanMd5'));
else
  % MATLAB
  docNode = com.mathworks.xml.XMLUtils.createDocument('PTTxPlanMd5');
end

  docNode.appendChild(docNode.createComment('TREATMENT PLAN CREATED WITH MATRAD'));

  PTTxPlanMd5 = docNode.getDocumentElement;
  PTTxPlanMd5.setAttribute('md5','noMD5');
  PTTxPlanMd5.setAttribute('xmlns:xsi','http://www.w3.org/2001/XMLSchema-instance');
  PTTxPlanMd5.setAttribute('xsi:noNamespaceSchemaLocation','RTT-PT-Plan.xsd');

  PTTxPlan = docNode.createElement('PTTxPlan');
  PTTxPlanMd5.appendChild(PTTxPlan);

  Beam = docNode.createElement('Beam');
  Beam.setAttribute('uid','bee035c5-03f6-4e9c-94b9-31f0fc484db1');
  PTTxPlan.appendChild(Beam);

  RstFormat = docNode.createElement('RstFormat');
  RstFormat.appendChild(docNode.createTextNode('PT_2004'));
  Beam.appendChild(RstFormat);

  Patient = docNode.createElement('Patient');
  Patient.setAttribute('id','PT-2004-01');
  Patient.setAttribute('name','unknown');
  Patient.setAttribute('sex','unknown');
  Patient.setAttribute('birthDate','2001-01-01');
  Beam.appendChild(Patient);

  TxInitiation = docNode.createElement('TxInitiation');
  TxInitiation.setAttribute('therapist','None');
  TxInitiation.setAttribute('dateTime','2007-01-23T13:52:27.2343750+01:00');
  Beam.appendChild(TxInitiation);

  TxRoom = docNode.createElement('TxRoom');
  TxRoom.setAttribute('name','Room1Fixed90');
  if strcmp(XMLplan(beamNb).ion,'1H')
    TxRoom.setAttribute('projectile','PROTON');
    TxRoom.setAttribute('charge','');
    TxRoom.setAttribute('mass','');
    TxRoom.setAttribute('atomicNumber','');

  elseif strcmp(XMLplan(beamNb).ion,'12C')
    TxRoom.setAttribute('projectile','ION');
    TxRoom.setAttribute('charge','6');
    TxRoom.setAttribute('mass','12');
    TxRoom.setAttribute('atomicNumber','6');
  end
  Beam.appendChild(TxRoom);

  BAMS = docNode.createElement('BAMS');
  if strcmp(XMLplan(beamNb).ion,'1H')
    BAMS.setAttribute('rippleFilter','254'); % value obtained from a xml plan sample from HIT
    BAMS.setAttribute('rangeShifter','254'); % value obtained from a xml plan sample from HIT
  elseif strcmp(XMLplan(beamNb).ion,'12C')
    BAMS.setAttribute('rippleFilter','3');   % value obtained from a xml plan sample from HIT
    BAMS.setAttribute('rangeShifter','254'); % value obtained from a xml plan sample from HIT
  end
  BAMS.setAttribute('rangeShifterDistance','0'); % value obtained from a xml plan sample from HIT
  Beam.appendChild(BAMS);

  TxTable = docNode.createElement('TxTable');
  TxTable.setAttribute('roll','0');         % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('pitch','0');        % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('lateral','0');      % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('longitudinal','0'); % value obtained from a xml plan sample from HIT
  aValue = sprintf('%d',XMLplan(beamNb).couchAngle);
  TxTable.setAttribute('isocentricAngle',aValue);
  TxTable.setAttribute('vertical','0');     % value obtained from a xml plan sample from HIT
  Beam.appendChild(TxTable);

  Gantry = docNode.createElement('Gantry');
  aValue = sprintf('%d',XMLplan(beamNb).gantryAngle);
  Gantry.setAttribute('angle',aValue);
  Beam.appendChild(Gantry);

  for iesNb = 1:length(XMLplan(beamNb).IES)

    spotNb = 1;
    
    IES = docNode.createElement('IES');
    aValue = sprintf('%d',iesNb);
    IES.setAttribute('number',aValue); 

    aValue = sprintf('%.2f',XMLplan(beamNb).IES(iesNb).energy);
    IES.setAttribute('energy',aValue);

    aValue = sprintf('%.1f',XMLplan(beamNb).IES(iesNb).focus);
    IES.setAttribute('focus',aValue);

    voxelMatrix = sortrows(XMLplan(beamNb).IES(iesNb).voxelMatrix,[1,2]);

      while length(voxelMatrix)

          Voxel = docNode.createElement('Voxel');
          aValue = sprintf('%f',voxelMatrix(1,1));
          Voxel.setAttribute('x',aValue);
          aValue = sprintf('%f',voxelMatrix(1,2));
          Voxel.setAttribute('y',aValue);
          aValue = sprintf('%d',voxelMatrix(1,3));
          Voxel.setAttribute('particles',aValue);
          IES.appendChild(Voxel);

          XMLplan(beamNb).IES(iesNb).rasterScanPath.x(spotNb,1) = voxelMatrix(1,1);
          XMLplan(beamNb).IES(iesNb).rasterScanPath.y(spotNb,1) = voxelMatrix(1,2);
          XMLplan(beamNb).IES(iesNb).rasterScanPath.particles(spotNb,1) = voxelMatrix(1,3);

          spotNb = spotNb + 1;
          
          tmpY = voxelMatrix(1,2);
          voxelMatrix(1,:) = [];

          if length(voxelMatrix)
            if (voxelMatrix(1,2) < tmpY)
              voxelMatrix = sortrows(voxelMatrix,[1,-2]);
            else
              voxelMatrix = sortrows(voxelMatrix,[1,2]);
            end
          end
      end
      Beam.appendChild(IES);
  end

  xmlwrite(filename,docNode);

end % loop over beams
