function resultGUI = findUnderdoseVolume(cst,ct,resultGUI)

% find PTV
for i = 1:size(cst,1)
    if strcmp(cst{i,2},'PTV') || strcmp(cst{i,2},'ITV')
        ptv = cst{i,4}{1};
    end
end

% find underdose in PTV for homogeneous lung
underIxHomo = find(resultGUI.matRadRecalc_RBExDose(ptv) < 2*.95);
underCubeHomo = zeros(ct.cubeDim);
underCubeHomo(ptv(underIxHomo)) = 1;
resultGUI.underHomo = underCubeHomo;

% find underdose in PTV for heterogeneous lung
underIxHetero = find(resultGUI.matRadHetero_RBExDose(ptv) < 2*.95);
underCubeHetero = zeros(ct.cubeDim);
underCubeHetero(ptv(underIxHetero)) = 1;
resultGUI.underHetero = underCubeHetero;

% find underdose in PTV for homogeneous lung only
underIxHomoOnly = setdiff(underIxHomo,underIxHetero);
underCubeHomoOnly = zeros(ct.cubeDim);
underCubeHomoOnly(ptv(underIxHomoOnly)) = 1;
resultGUI.underHomoOnly = underCubeHomoOnly;

% find underdose in PTV for heterogeneous lung only
underIxHeteroOnly = setdiff(underIxHetero,underIxHomo);
underCubeHeteroOnly = zeros(ct.cubeDim);
underCubeHeteroOnly(ptv(underIxHeteroOnly)) = 1;
resultGUI.underHeteroOnly = underCubeHeteroOnly;

