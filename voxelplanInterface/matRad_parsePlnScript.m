function pln = matRad_parsePlnScript(filename)

h = fopen(filename,'r');

while ~feof(h)
   
    currLine = fgetl(h);
    
    if strncmp(currLine,'pln.SAD',7)
        pln.SAD = sscanf(currLine,'pln.SAD = %f;');
    elseif strncmp(currLine,'pln.bixelWidth',14)
        pln.bixelWidth = sscanf(currLine,'pln.bixelWidth = %f;');
    elseif strncmp(currLine,'pln.gantryAngles',16)
        pln.gantryAngles = str2num(currLine(find(currLine=='[',1)+1:find(currLine==']',1)-1));
    elseif strncmp(currLine,'pln.couchAngles',15)
        pln.couchAngles = str2num(currLine(find(currLine=='[',1)+1:find(currLine==']',1)-1));
    elseif strncmp(currLine,'pln.radiationMode',17)
        pln.radiationMode = sscanf(currLine,'pln.radiationMode = ''%[^'']');
    elseif strncmp(currLine,'pln.bioOptimization',19)
        pln.bioOptimization = sscanf(currLine,'pln.bioOptimization = ''%[^'']');
    elseif strncmp(currLine,'pln.numOfFractions',18)
        pln.numOfFractions = sscanf(currLine,'pln.numOfFractions = %f;');
    elseif strncmp(currLine,'pln.runSequencing',17)
        pln.runSequencing = str2num(sscanf(currLine,'pln.runSequencing = %[^;]'));
    elseif strncmp(currLine,'pln.runDAO',10)
        pln.runDAO = str2num(sscanf(currLine,'pln.runDAO = %[^;]'));
    end   

end

fclose(h);
