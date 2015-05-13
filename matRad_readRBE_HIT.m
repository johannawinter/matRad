clc
clear
close all

Spectra = {'hydrogen','helium','lithium','beryllium','bor','carbon','nitrogen','oxygen','fluor','neon'};

path = 'E:\TRiP98DATA_HIT-20131120\RBE';
folderInfo = dir(path);
CntFiles = 1;

for i = 1:length(folderInfo)
   
    if folderInfo(i).isdir ~=1 && folderInfo(i).bytes>0
       
        fileID = fopen([path filesep folderInfo(i).name]);
        data = textscan(fileID,'%s');
        data = data{1,1};
      
        FlagParse = true;
        Cnt = 1;
        % parse header
        while FlagParse
           
           if ~isempty(strfind(data{Cnt,1},'projectile'))
               FlagParse = false;
               break;
           elseif ~isempty(strfind(data{Cnt,1},'mean'))
               meta(CntFiles).(strrep(data{Cnt,1},'!','')) = '';
                 Cnt = Cnt +1;
           else
               % parse input as number or as string
               if isempty(str2num(data{Cnt+1,1}))
                   val = data{Cnt+1,1};
               else
                   val = str2num(data{Cnt+1,1});
               end
               meta(CntFiles).(strrep(data{Cnt,1},'!','')) = val;
               Cnt = Cnt +2;
           end

        end
      
        % parse projectile headers
        for SpecCnt = 1:length(Spectra)
            %meta(CntFiles).(strrep(data{Cnt,1},'!','')){1,1} = data{Cnt+1,1
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.projectile = data{Cnt+1,1};
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.A = SpecCnt*2;
            meta(CntFiles).(Spectra{SpecCnt}){1,1}.Z = SpecCnt;
            isNumber = true;
            CntDat = Cnt + 6;
            CntStruct = 1;
            while isNumber
                 if CntDat>length(data)
                    break;
                 end
                 if isempty(str2num(data{CntDat,1})) || isempty(str2num(data{CntDat+1,1})) 
                     isNumber = false;
                     break;
                 else
                     s(CntStruct).Energy = str2num(data{CntDat,1});
                     s(CntStruct).RBE = str2num(data{CntDat+1,1});
                     CntStruct = CntStruct+1;
                     CntDat = CntDat+2;
                 end
            end

            meta(CntFiles).(Spectra{SpecCnt}){1,2} = s;
            Cnt = CntDat;
        end
        
     CntFiles = CntFiles+1;
    end
end


%% plot RBE spectra of specific cell type;
for j = 1:23
    CellType = j;
    figure,

    for i = 1:length(Spectra)
        vX = [meta(CellType).(Spectra{i}){1,2}(:).Energy];
        vY = [meta(CellType).(Spectra{i}){1,2}(:).RBE];
        plot(vX,vY,'Linewidth',3),hold on
    end

    str = sprintf('celltype: alpha_x: %f and beta_x: %f',meta(CellType).alpha,meta(CellType).beta);
    legend(Spectra),xlabel('energy [MeV/u]'),ylabel('RBE'),grid on
    title(str);
    set(gca,'FontSize',16)
end


%%
% clc
% clear
% close all
% 
% pathSpec = 'E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV08000.spc';
% 
% fileID = fopen(pathSpec);
% data = fread(fileID);
% 
% tgspcread


