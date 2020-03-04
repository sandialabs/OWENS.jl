function processAeroLoads()
    d1 = dir('*.geom');

    d2 = dir('*_ElementData.csv');

    d3 = dir('*.bld');

    d4 = dir('*.el');

    d5 = dir('*.ort');

    d6 = dir('*.mesh');

    [timeArray,ForceValHist,ForceDof] = mapCactusLoadsFile(d1.name,d2.name,d3.name,d4.name,d5.name,d6.name);

    clear d1 d2 d3 d4 d5 d6

    save aeroLoads
end