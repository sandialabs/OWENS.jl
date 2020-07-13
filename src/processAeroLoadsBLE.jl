function [aeroLoads] = processAeroLoadsBLE(aeroLoadsFile, OWENSfile)

    aeroLoadsFile_root = aeroLoadsFile(1:end-16); #cut off the _ElementData.csv
    OWENSfile_root = OWENSfile(1:end-6); #cut off the .owens

    d1 = [aeroLoadsFile_root '.geom'];
    d2 = [aeroLoadsFile_root '_ElementData.csv'];
    d3 = [OWENSfile_root '.bld'];
    d4 = [OWENSfile_root '.el'];
    d5 = [OWENSfile_root '.ort'];
    d6 = [OWENSfile_root '.mesh'];

    [timeArray,ForceValHist,ForceDof] = mapCactusLoadsFile(d1,d2,d3,d4,d5,d6);
    aeroLoads.timeArray = timeArray;
    aeroLoads.ForceValHist = ForceValHist;
    aeroLoads.ForceDof = ForceDof;

    # method constrained by not passing the .mat filename
    # save('aeroLoads.mat','timeArray','ForceValHist','ForceDof');
    # disp('New aeroLoads.mat file saved.')

end
