function processAeroLoadsBLE(CACTUSfileRoot, OWENSfileRoot, outputAeroFileName)

if 1
    d1 = dir([CACTUSfileRoot '.geom']);
    d2 = dir([CACTUSfileRoot '_ElementData.csv']);
    d3 = dir([OWENSfileRoot '.bld']);
    d4 = dir([OWENSfileRoot '.el']);
    d5 = dir([OWENSfileRoot '.ort']);
    d6 = dir([OWENSfileRoot '.mesh']);
else
    d1 = dir('*.geom');
    d2 = dir('*_ElementData.csv');
    d3 = dir('*.bld');
    d4 = dir('*.el');
    d5 = dir('*.ort');
    d6 = dir('*.mesh');
end

[timeArray,ForceValHist,ForceDof] = mapCactusLoadsFile(d1.name,d2.name,d3.name,d4.name,d5.name,d6.name);

clear d1 d2 d3 d4 d5 d6
clear CACTUSfileRoot OWENSfileRoot


if 0 % preferred method
    save TEMP_aeroLoads
    movefile('TEMP_aeroLoads.mat',[outputAeroFileName '.mat'])
else % method constrained by not passing the .mat filename
    save aeroLoads
    disp('New aeroLoads.mat file saved.')
end

end