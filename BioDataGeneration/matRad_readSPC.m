function SPC = matRad_calcParticleDose(path)

% TRiP98 spectra file format
% Spectra (SPC) files are binary files containing energy spectra and related histograms of 
% the various particles created when an ion undergoes nuclear interactions with the traversed matter. 
% There is one file per initial beam energy. Spectra are stored as a function of depth in water. 
% The file is organized in a similar fashion as a TIFF file, i.e. each data entry is preceeded by a so-called 
% tag field which identifies the following data item by a unique number and contains the length of that item. 
% This way a reader progam can skip irrelevant or unknown data items. With a few exceptions the 
% data items (including the tags) are binary, and thus can be directly read only on a CPU architecture 
% with the same byte order. This was implemented deliberately, since spectra are meant to be memory-mapped 
% read-only into CPU address space, which disallows byte swap operations. Memory mapping also implies that 
% structured data must be written padded to 8-byte boundaries, otherwise data access time penalties and/or
% access violations may occur. Hence 4-byte integers are expanded to 8 bytes and character strings always end 
% on 8-byte boundaries, with binary zeroes padded. Floating point numbers in TRiP98 are always 8-byte double 
% precision IEEE numbers anyway.
% 
% A tag entry is defined as
% 
% struct STRPSPCBTAG { unsigned long ulTag; unsigned long ulLen; };
% 
% which is a sequence of two 4-byte integers. The Len field specifies the length (in bytes) of the data 
%item to follow, excluding the tag entry. The currently supported Tag codes are:
% 
% enum {
%        TRPSPCBTAG_FILETYPE   =1,      /* header info */
%        TRPSPCBTAG_FILEVERSION=2,
%        TRPSPCBTAG_FILEDATE   =3,
%        TRPSPCBTAG_TARGNAME   =4,
%        TRPSPCBTAG_PROJNAME   =5,
%        TRPSPCBTAG_B          =6,
%        TRPSPCBTAG_P          =7,
%        TRPSPCBTAG_N          =8,
%        TRPSPCBTAG_NZ         =9,
%        TRPSPCDTAG_Z          =10,
%        TRPSPCDTAG_N          =11,       
%        TRPSPCDTAG_NS         =12,
%        TRPSPCDTAG_S          =13,
%        TRPSPCDTAG_CUM        =14,
%        TRPSPCDTAG_NC         =15,
%        TRPSPCDTAG_NE         =16,
%        TRPSPCDTAG_E          =17,
%        TRPSPCDTAG_EREF       =18,
%        TRPSPCDTAG_HISTO      =19,
%        TRPSPCDTAG_RUNNINGSUM =20,
% };
% 
% The file format starts as follows:
% 
% BTAG_FILETYPE    <filetype> 
% BTAG_FILEVERSION <fileversion> 
% BTAG_FILEDATE    <filedate> 
% BTAG_TARGNAME    <targname>
% BTAG_PROJNAME    <projname>
% BTAG_B           <B>        # double; beam energy [MeV/u]
% BTAG_P           <P>        # double; peak position [g/cm**2]
% BTAG_N           <N>        # double; normalization, usually =1
% BTAG_NZ          <nz>       # 8-byte unsigned integer; number of depth steps.
% <nz> depth data blocks
% 
% where <filetype> is an 80-byte ASCII character string starting with "SPCM" or "SPCI", 
%specifying big ("Motorola") or little ("Intel") endian byte order, respectively.
% <fileversion> is an 80-byte ASCII character string specifying the file format version 
%as yyyymmdd. 19980704 is the "classic" format (fixed energy range), whereas 20020107 is 
%reserved for future possible variable energy range.
% <filedate> is an 80-byte ASCII character string with the file creation date as returned 
%by the ctime() function. (<dow> <mmm> <dd> <hh>:<mm>:<ss> <yyyy>)
% <targname> and <projname> are the names of target ("H2O") and projectile ("12C6"), respectively.
%Since both can have any length, they are padded to the right with binary zeroes up to the next 8-byte boundary.
% A depth data block is organized as follows:
% 
% DTAG_Z <P>        # double; depth [g/cm**2]
% DTAG_N <N>        # double; normalization for this depth step, usually =1.
% DTAG_NS <nS>      # 8-byte unsigned integer; number of particle species.
% <nS> species data blocks
% 
% A species data block is organized as follows:
% 
% DTAG_S    <ZA>     # double Z, double A, long Z, long A;
% DTAG_CUM  <Cum>    # double
% DTAG_NC   <nC>     # 8-byte unsigned integer; 
% DTAG_NE   <nE>     # 8-byte unsigned integer; number of energy bins for this species
% {
% DTAG_EREF <lSRef>  # 8-byte unsigned integer; 
% |
% DTAG_E <E[nE+1]>   # double; energy bin values 
% }
% DTAG_HISTO <H[nE]>           # double; spectrum bin values 
% DTAG_RUNNINGSUM <Cum[nE+1]>  # double; running cumulated spectrum bin values 
% 
% The scalar Cum value is the cumulated number (running sum) of fragments, i.e. the species sum over Cum[nE]. 
%This number may exceed 1, since for an incoming primary particle many secondaries are created.
% <nC> is reserved for later use, so that lateral scattering for each fragment can be included. At present nC=0.
% If a species data block is flagged as EREFCOPY, its energy bins are not stored, but rather a reference
%index <lSRef> (0..<nS>-1) to a species with identical energy binning. This helps to reduce space requirements, 
%since fragment spectra may share the same energy bins. If the energy bins for the species under consideration are 
%not a reference copy, <nE>+1 double precision bin border values are stored.
% The <H[]> are the spectrum contents divided by the bin width.
% The <Cum[]> are the running integrated values of <H[]>
% 
% The usual file name extension is .spc.
% The usual naming convention is <pp>.<tt>.<uuu><eeeee>.spc, where <pp> denotes the projectile, 
%<tt> the target material, <uuu> the unit ( keV, MeV, GeV) and <eeeee> the energy in these units, 
%with the decimal point after the middle digit. Example: 12C.H2O.MeV27000.spc refers to 270 MeV/u. 

% information extracted from:
% http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html

if isempty(path)
    path = ('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV09000.spc');
end


h = waitbar(0,'Initializing waitbar...');

%% create Tag lookup table
TagMAP = containers.Map;
TagMAP('1')  = 'filetype';
TagMAP('2')  = 'fileversion';
TagMAP('3')  = 'filedate';
TagMAP('4')  = 'target';
TagMAP('5')  = 'projectile';
TagMAP('6')  = 'energy';
TagMAP('7')  = 'peakPos';
TagMAP('8')  = 'normalization';
TagMAP('9')  = 'numDepthSteps';
TagMAP('10') = 'depths';
TagMAP('11') = 'depthNormalization';
TagMAP('12') = 'numParticles';
TagMAP('13') = 'typeOfParticle';
TagMAP('14') = 'CumSumFragements';
TagMAP('15') = 'nC';
TagMAP('16') = 'numOfEnergyBin';
TagMAP('17') = 'energyBinValue';
TagMAP('18') = 'Eref';
TagMAP('19') = 'Hist';
TagMAP('20') = 'RunningSum';

sParticleNames = {'H','He','Li','Be','B','C'};

SPC = struct;
% second argument 'l' or 'b' determiens endian - little or big
fid = fopen(path,'rb','l');

%% parse meta data
LengthMetaData = 9;
for i=1:LengthMetaData

    [Tag, TagLength] = findNextTag(fid);
          
      switch Tag
        % first 5 entries are characters  
        case {1,2,3,4,5}
            SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'char');    
        % read double
        case {6,7,8}
            SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'double');          
        % read 8-byte unsigned integer
        case  9
            SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'integer'); 
       end
         
end



%% parse spectrum 

ReferenceBlockIdx = 0;
depthStep = 0;
CurrPos = ftell(fid);
fseek(fid,0,'eof');
EndPos = ftell(fid);
fseek(fid,CurrPos,'bof');

while true
    
   [Tag, TagLength] = findNextTag(fid,EndPos);
     
        switch Tag
            
            % read depth as double - Tag 10 indicates new depth block
            case 10
                 depthStep = depthStep + 1;
                 Process = double(depthStep)/double(SPC.(TagMAP('9')));
                 waitbar(Process,h,'parsing depth blocks');
                 
                 SPC.Data(depthStep).depthStep     = depthStep;
                 SPC.Data(depthStep).(TagMAP('4')) = SPC.(TagMAP('4'));
                 SPC.Data(depthStep).(TagMAP('5')) = SPC.(TagMAP('5'));
                 SPC.Data(depthStep).(TagMAP('6')) = SPC.(TagMAP('6'));
                 SPC.Data(depthStep).(TagMAP('7')) = SPC.(TagMAP('7'));
                 SPC.Data(depthStep).(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'double'); 
            
            % read normalization as double
            case {11,14}   
                 SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'double');

            % read number of particle species as integer
            case  12
                SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'integer');
                
            % read atomic number z and mass number a as double and long 
            case 13
                
                % Tag 13 indicates a new depth block
                if isfield(SPC,TagMAP('13')) 
                    if size(SPC.(TagMAP(num2str(Tag))),2) < SPC.(TagMAP('12'))
                        SPC.(TagMAP('13')) = horzcat(SPC.(TagMAP('13')),readSpcValue(fid,16,'double'));
                    else
                        % read it out but dont store it
                        tmp = readSpcValue(fid,16,'double');
                    end
                else
                    SPC.(TagMAP('13')) = readSpcValue(fid,16,'double');
                end
              
                A_Z = readSpcValue(fid,8,'long');
                
                if sum(A_Z == [1;2])==2
                    CurrentParticle = 'H';
                elseif sum(A_Z == [2;4])==2
                    CurrentParticle = 'He';
                elseif sum(A_Z == [3;6])==2
                    CurrentParticle = 'Li';
                elseif sum(A_Z == [4;8])==2
                    CurrentParticle = 'Be';
                elseif sum(A_Z == [5;10])==2
                    CurrentParticle = 'B';
                elseif sum(A_Z == [6;12])==2
                    CurrentParticle = 'C';
                end
                
            case {15,16}
                 % number of energy bins for a specific species
                SPC.(TagMAP(num2str(Tag))) = readSpcValue(fid,TagLength,'long');
                
            case 17
                % parse low and high energy borders of each bin and
                % calculated mid energy
                NumEnergies =  SPC.(TagMAP('16'));
                Elow = zeros(1,NumEnergies);
                Ehigh = zeros(1,NumEnergies);
                Emid = zeros(1,NumEnergies);
                dE = zeros(1,NumEnergies);
                
                for idx = 1:NumEnergies + 1
                    if idx == 1
                       Elow(idx)  = readSpcValue(fid,TagLength,'double');
                    else
                       Elow(idx)    = readSpcValue(fid,TagLength,'double');
                       Ehigh(idx-1) = Elow(idx);
                       Emid(idx-1)  =  sqrt(Ehigh(idx-1) * Elow(idx-1));
                       %Emid  = (Ehigh(idx-1) + Elow(idx-1))/2;
                       dE(idx-1) = Ehigh(idx-1) - Elow(idx-1);
                    end
                end
                
                Elow(end) = [];
                SPC.Data(depthStep).(CurrentParticle).Elow = Elow;
                SPC.Data(depthStep).(CurrentParticle).Emid = Emid;
                SPC.Data(depthStep).(CurrentParticle).Ehigh = Ehigh;
                SPC.Data(depthStep).(CurrentParticle).dE = dE;
                SPC.Data(depthStep).(CurrentParticle).ReferenceBlockIdx = ReferenceBlockIdx;
                ReferenceBlockIdx = ReferenceBlockIdx + 1;
               
            case 18
                
               % energy bins are reused to save some space
               % reference to a species with identical energy binning is
               % provided  (0..<nS>-1)

               RefIdx = readSpcValue(fid,8,'double');
               
               if RefIdx > ReferenceBlockIdx
                disp('reference to non existing energy block');
               end
               
               %resolve index
               DepthStepRef  = (RefIdx/SPC.(TagMAP('12'))) + 1;
               ParticleNoRef = (mod(RefIdx,SPC.(TagMAP('12')))) + 1; 
             
               % find depth block index
               SPC.Data(depthStep).(CurrentParticle).Elow = SPC.Data(DepthStepRef).(sParticleNames{ParticleNoRef}).Elow;
               SPC.Data(depthStep).(CurrentParticle).Emid = SPC.Data(DepthStepRef).(sParticleNames{ParticleNoRef}).Emid;
               SPC.Data(depthStep).(CurrentParticle).Ehigh = SPC.Data(DepthStepRef).(sParticleNames{ParticleNoRef}).Ehigh;
               SPC.Data(depthStep).(CurrentParticle).dE = SPC.Data(DepthStepRef).(sParticleNames{ParticleNoRef}).dE;
               
               SPC.Data(depthStep).(CurrentParticle).ReferenceBlockIdx = ReferenceBlockIdx;
               ReferenceBlockIdx = ReferenceBlockIdx + 1;
                
            case 19
                               
                NumEnergies = SPC.(TagMAP('16'));
                dNdE = readSpcValue(fid,8*NumEnergies,'double');
                SPC.Data(depthStep).(CurrentParticle).dNdE = dNdE';
                     
            case 20
                
                % read the number of particles per primary for each bin
                NumEnergies =  SPC.(TagMAP('16'))+1;
               
                for idx = 1:NumEnergies+1                 
                    if idx == 1
                      N(idx)  = readSpcValue(fid,8,'double');
                    else
                      N(idx) = readSpcValue(fid,8,'double') - sum([N(1:idx-1)]);   
                    end

                end
                
               SPC.Data(depthStep).(CurrentParticle).N = N(2:end);
        
            case -1
                % reached end of file
                close(h);
                break;
                
        end

end



end




function value = readSpcValue(fid,NumOfBytesToRead,ValueType)

    switch ValueType
        
        case 'double'
             value = uint8(fread(fid,NumOfBytesToRead));
             value = typecast(value,'double');

        case 'integer'
            value = uint8(fread(fid,NumOfBytesToRead));
            value(value==0)=[];
            
            if length(value) == 2
               value = value(1);
               %disp('dirty!')
            end
            
        case 'char'
            value = (fread(fid,NumOfBytesToRead));
            value(value==0)=[];
            value =(char(value))';
            
        case 'long'
             value = uint8(fread(fid,NumOfBytesToRead));
             value = typecast(value,'uint32');
             value(value==0)=[];
    end
end


function [Tag, Length] = findNextTag(fid,EndPos)

      while true
            
            Tag = fread(fid,1);
            if Tag > 0 
                break;
            end
            
            if ftell(fid)>=EndPos
               Tag = -1;
               break;
            end
      end
       
      Length = readSpcValue(fid,7,'integer');
      
end




