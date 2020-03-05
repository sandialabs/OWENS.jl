function processAeroLoadsBLE(CACTUSfileRoot, OWENSfileRoot, outputAeroFileName)

if 1
    d1 = [CACTUSfileRoot '.geom'];
    d2 = [CACTUSfileRoot '_ElementData.csv'];
    d3 = [OWENSfileRoot '.bld'];
    d4 = [OWENSfileRoot '.el'];
    d5 = [OWENSfileRoot '.ort'];
    d6 = [OWENSfileRoot '.mesh'];
else
    d1 = dir('*.geom');
    d2 = dir('*_ElementData.csv');
    d3 = dir('*.bld');
    d4 = dir('*.el');
    d5 = dir('*.ort');
    d6 = dir('*.mesh');
end

[timeArray,ForceValHist,ForceDof] = mapCactusLoadsFile(d1,d2,d3,d4,d5,d6);

%clear d1 d2 d3 d4 d5 d6
%clear CACTUSfileRoot OWENSfileRoot


if 0 % preferred method
    save TEMP_aeroLoads
    movefile('TEMP_aeroLoads.mat',[outputAeroFileName '.mat'])
else % method constrained by not passing the .mat filename
    save aeroLoads
    disp('New aeroLoads.mat file saved.')
end

end
