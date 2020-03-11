function processAeroLoadsBLE(CACTUSfileRoot, OWENSfileRoot)

    d1 = [CACTUSfileRoot '.geom'];
    d2 = [CACTUSfileRoot '_ElementData.csv'];
    d3 = [OWENSfileRoot '.bld'];
    d4 = [OWENSfileRoot '.el'];
    d5 = [OWENSfileRoot '.ort'];
    d6 = [OWENSfileRoot '.mesh'];

    [timeArray,ForceValHist,ForceDof] = mapCactusLoadsFile(d1,d2,d3,d4,d5,d6);

    % method constrained by not passing the .mat filename
    save('aeroLoads.mat','timeArray','ForceValHist','ForceDof');
    disp('New aeroLoads.mat file saved.')

end
