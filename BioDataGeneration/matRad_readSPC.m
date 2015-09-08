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

clc
clear all
close all
path = ('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV09000.spc');
% http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html

filename = path;		% hypothetical file

TagMAP = containers.Map;
TagMAP('1') = 'FILETYPE';
TagMAP('2') = 'FILEVERSION';
TagMAP('3') = 'FILEDATE';
TagMAP('4') = 'TARGET';
TagMAP('5') = 'PROJECTILE';
TagMAP('6') = 'BEAM_ENERGY_MeV_u';
TagMAP('7') = 'PEAK_POS_g_cm_2';
TagMAP('8') = 'NORMALIZATION';
TagMAP('9') = 'NUM_DEPTH';

spc = struct;


CurrPos = 0;
fid = fopen(filename,'rt','b');

for i=1:9
    % find unique tag number
    while true
        Tag = fread(fid,CurrPos+1);
        if Tag >0
            break;
        end
    end
    
    % find byte length
    while true
        TagLength = fread(fid,CurrPos+1);
        if TagLength >0
            break;
        end
    end
  % first 5 entries are characters   
  if i<6
    Tmp = fread(fid,TagLength);
    Tmp(Tmp==0) = [];
    spc.(TagMAP(num2str(Tag))) = (char(Tmp))';
    
  else

    switch Tag
        % read double
        case {6,7,8,10}
            fseek(fid, 3, 0);
            Tmp = fread(fid,TagLength);
            spc.(TagMAP(num2str(Tag))) = typecast(uint8(Tmp),'double');
            
        % read 8-byte unsigned integer
        case  9
            Tmp = fread(fid,TagLength);
            Tmp(Tmp==0) = [];
            spc.(TagMAP(num2str(Tag))) = Tmp;
    end
            
    
  end
end


%% start parsing depth data blocks
    %% parse each depth block

% create
spc.data(1).depthStep = 1;
spc.data(1).H.Elow = 1;
spc.data(1).H.Emid =  1;
spc.data(1).H.Ehigh =  1;
spc.data(1).H.dE = 1;
spc.data(1).H.dNdE = 1;
spc.data(1).H.N = 1;

% spc.data(1).He = 
% spc.data(1).Li =
% spc.data(1).Be =
% spc.data(1).B =
% spc.data(1).C =
