function  matRad_XlsSpectra2Mat(pathSPC,pathSaving)
%% read spc file
Directory = dir([pathSPC filesep '*.xlsx']);

sParticles = {'H','He','Li','Be','B','C'};
sParticlesNo = {'1002','2004','3006','4008','5010','6012'};


for CntFile = 1:length(Directory)

    [~,~,raw] = xlsread([pathSPC filesep Directory(CntFile).name]);
    % skip the first column and first row
    raw = raw(2:end,2:end);
    % convert the data to numbers - takes quite a bit
    rawNum = str2double(raw(:,1:10));
    Cnt = 1;
    for iDepth=1:79
        % save meta information
        SPC(iDepth,1).depthStep = iDepth;
        SPC(iDepth,1).depths = rawNum(Cnt,2);
        SPC(iDepth,1).projectile = '12C';
        SPC(iDepth,1).target = 'H2O';
        SPC(iDepth,1).energy = raw{2,13};
        SPC(iDepth,1).peakPos = raw{2,14};

        % extract the spectra data for current depth step
        for IdxPart = 1:length(sParticles)
           InnerCnt = 1;
           
            while true

                if Cnt > length(raw)
                    break;
                end

               if strcmp(sParticlesNo(IdxPart),raw{Cnt,3})
                    Elow{InnerCnt} = rawNum(Cnt,4);
                    Emid{InnerCnt} = rawNum(Cnt,5);
                    Ehigh{InnerCnt} = rawNum(Cnt,6);
                    dE{InnerCnt} = rawNum(Cnt,7);
                    dNdE{InnerCnt} = rawNum(Cnt,8);
                    N{InnerCnt} = rawNum(Cnt,9);
                    InnerCnt = InnerCnt+1;
                    Cnt = Cnt +1;
               else
                   SPC(iDepth,1).(sParticles{IdxPart}).Elow = cell2mat(Elow);
                   Elow = [];
                   SPC(iDepth,1).(sParticles{IdxPart}).Emid = cell2mat(Emid);
                   Emid = [];
                   SPC(iDepth,1).(sParticles{IdxPart}).Ehigh = cell2mat(Ehigh);
                   Ehigh=[];
                   SPC(iDepth,1).(sParticles{IdxPart}).dE = cell2mat(dE);
                   dE=[];
                   SPC(iDepth,1).(sParticles{IdxPart}).dNdE = cell2mat(dNdE);
                   dNdE=[];
                   SPC(iDepth,1).(sParticles{IdxPart}).N = cell2mat(N);
                   N=[];
                   % stop while loop
                   break
               end 
            end
        end 
    end
    
    
    % set last values to zero
    SPC(79).C.Elow = 0;
    SPC(79).C.Emid =0;
    SPC(79).C.Ehigh =0;
    SPC(79).C.dE =0;
    SPC(79).C.dNdE =0;
    SPC(79).C.N =0;
    
    if sum(strcmp(SPC(1).energy,{'80','90'}))>0
        VarName = ['C12spc' '0' SPC(1).energy];
    else
        VarName = ['C12spc' num2str(SPC(1).energy)];
    end
    
    eval([VarName ' = ' 'SPC']);
    save([pathSaving filesep VarName],VarName);

end

