clc
clear all
close all
path = ('E:\TRiP98DATA_HIT-20131120\SPC\12C\RF3MM\FLUKA_NEW3_12C.H2O.MeV09000.spc');


filename = path;		% hypothetical file
segsize = 100000;

%fid = fopen(filename,'r','b');


Encoding = {'Big5','ISO-8859-1','windows-932','GB2312','ISO-8859-2',...
    'windows-936','EUC-KR','ISO-8859-3','windows-949','EUC-JP','ISO-8859-4',...
	'windows-950','GBK','ISO-8859-9','windows-1250','KSC_5601','ISO-8859-13','windows-1251',...
    'Macintosh','ISO-8859-15','windows-1252','Shift_JIS','windows-1253',...
    'US-ASCII','windows-1254','UTF-8','windows-1257'};

fid = fopen(filename,'r','l','UTF-8');
StartPos = fseek(fid,0,'bof');
EndPos = fseek(fid, 0, 'eof');
tags = 1;
Maxtags = 20000;
mTag = zeros(4,Maxtags);

FirstLine = fgetl(fid)
SecondLine = fgetl(fid)
ThirdLine = fgetl(fid)

test = fread(fid)
    
fclose(fid);
for i = 1:length(Encoding)
    
    Byte = fread(fid);
    for j = 1:5000
     testi{j}=char(Byte(j));
    end
end

while ~feof(fid)
    Byte = fread(fid);
    %fseek(fid, 2, 0);
    for i = 1:length(Encoding)
     unicodestr{i} = native2unicode(Byte(598),Encoding{1,i});
    end
   
end
    
fclose(fid);