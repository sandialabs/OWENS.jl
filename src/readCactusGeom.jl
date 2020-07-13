function [cactusGeom] = readCactusGeom(geomfn)

fid = fopen(geomfn,'r');

NBlade = int32(processLine(fid));
[cactusGeom.NBlade] = NBlade(1);
NStrut = int32(processLine(fid));
[cactusGeom.NStrut] = NStrut(1);
[cactusGeom.RotN] = processLine(fid);
[cactusGeom.RotP] = processLine(fid);
[cactusGeom.RefAR] = processLine(fid);
[cactusGeom.RefR] = processLine(fid);

myfgetl(fid); # skip a line

myfgetl(fid); # skip a line
NElem =  processLine(fid);
blade = struct('NElem',0.0,...
    'FlipN',0.0,...
    'QCx',zeros(NElem(1,1)+1,1),...
    'QCy',zeros(NElem(1,1)+1,1),...
    'QCz',zeros(NElem(1,1)+1,1),...
    'tx',zeros(NElem(1,1)+1,1),...
    'ty',zeros(NElem(1,1)+1,1),...
    'tz',zeros(NElem(1,1)+1,1),...
    'CtoR',zeros(NElem(1,1),1),...
    'PEx',zeros(NElem(1,1),1),...
    'PEy',zeros(NElem(1,1),1),...
    'PEz',zeros(NElem(1,1),1),...
    'tEx',zeros(NElem(1,1),1),...
    'tEy',zeros(NElem(1,1),1),...
    'tEz',zeros(NElem(1,1),1),...
    'nEx',zeros(NElem(1,1),1),...
    'nEy',zeros(NElem(1,1),1),...
    'nEz',zeros(NElem(1,1),1),...
    'sEx',zeros(NElem(1,1),1),...
    'sEy',zeros(NElem(1,1),1),...
    'sEz',zeros(NElem(1,1),1),...
    'ECtoR',zeros(NElem(1,1),1),...
    'EAreaR',zeros(NElem(1,1),1),...
    'iSect',zeros(NElem(1,1),1));

cactusGeom.blade = repmat(blade,1,NBlade(1));
for i=1:NBlade(1,1)

    if i ~= 1
        myfgetl(fid); # skip a line
        myfgetl(fid); # skip a line
    end

    cactusGeom.blade(1,i).NElem = NElem(1,1);
    FlipN = processLine(fid);
    cactusGeom.blade(1,i).FlipN(1,1) = FlipN(1,1);

    QCx =  processLine(fid);
    cactusGeom.blade(1,i).QCx(1:NElem(1,1)+1,1) = QCx(1:NElem(1,1)+1,1);
    QCy =  processLine(fid);
    cactusGeom.blade(1,i).QCy(1:NElem(1,1)+1,1) = QCy(1:NElem(1,1)+1,1);
    QCz =  processLine(fid);
    cactusGeom.blade(1,i).QCz(1:NElem(1,1)+1,1) = QCz(1:NElem(1,1)+1,1);

    tx =  processLine(fid);
    cactusGeom.blade(1,i).tx(1:NElem(1,1)+1,1) = tx(1:NElem(1,1)+1,1);
    ty =  processLine(fid);
    cactusGeom.blade(1,i).ty(1:NElem(1,1)+1,1) = ty(1:NElem(1,1)+1,1);
    tz =  processLine(fid);
    cactusGeom.blade(1,i).tz(1:NElem(1,1)+1,1) = tz(1:NElem(1,1)+1,1);

    CtoR =  processLine(fid);
    cactusGeom.blade(1,i).CtoR(1:NElem(1,1),1) = CtoR(1:NElem(1,1),1);

    PEx =  processLine(fid);
    cactusGeom.blade(1,i).PEx(1:NElem(1,1),1) = PEx(1:NElem(1,1),1);
    PEy =  processLine(fid);
    cactusGeom.blade(1,i).PEy(1:NElem(1,1),1) = PEy(1:NElem(1,1),1);
    PEz =  processLine(fid);
    cactusGeom.blade(1,i).PEz(1:NElem(1,1),1) = PEz(1:NElem(1,1),1);

    tEx =  processLine(fid);
    cactusGeom.blade(1,i).tEx(1:NElem(1,1),1) = tEx(1:NElem(1,1),1);
    tEy =  processLine(fid);
    cactusGeom.blade(1,i).tEy(1:NElem(1,1),1) = tEy(1:NElem(1,1),1);
    tEz =  processLine(fid);
    cactusGeom.blade(1,i).tEz(1:NElem(1,1),1) = tEz(1:NElem(1,1),1);

    nEx =  processLine(fid);
    cactusGeom.blade(1,i).nEx(1:NElem(1,1),1) = nEx(1:NElem(1,1),1);
    nEy =  processLine(fid);
    cactusGeom.blade(1,i).nEy(1:NElem(1,1),1) = nEy(1:NElem(1,1),1);
    nEz =  processLine(fid);
    cactusGeom.blade(1,i).nEz(1:NElem(1,1),1) = nEz(1:NElem(1,1),1);

    sEx =  processLine(fid);
    cactusGeom.blade(1,i).sEx(1:NElem(1,1),1) = sEx(1:NElem(1,1),1);
    sEy =  processLine(fid);
    cactusGeom.blade(1,i).sEy(1:NElem(1,1),1) = sEy(1:NElem(1,1),1);
    sEz =  processLine(fid);
    cactusGeom.blade(1,i).sEz(1:NElem(1,1),1) = sEz(1:NElem(1,1),1);

    ECtoR =  processLine(fid);
    cactusGeom.blade(1,i).ECtoR(1:NElem(1,1),1) = ECtoR(1:NElem(1,1),1);
    EAreaR =  processLine(fid);
    cactusGeom.blade(1,i).EAreaR(1:NElem(1,1),1) = EAreaR(1:NElem(1,1),1);
    iSect =  processLine(fid);
    cactusGeom.blade(1,i).iSect(1:NElem(1,1),1) = iSect(1:NElem(1,1),1);
end


myfgetl(fid); # skip line
NElem =  processLine(fid);
strut = struct('NElem',0.0,...
    'TtoC',0.0,...
    'MCx',zeros(NElem(1,1)+1,1),...
    'MCy',zeros(NElem(1,1)+1,1),...
    'MCz',zeros(NElem(1,1)+1,1),...
    'CtoR',zeros(NElem(1,1)+1,1),...
    'PEx',zeros(NElem(1,1),1),...
    'PEy',zeros(NElem(1,1),1),...
    'PEz',zeros(NElem(1,1),1),...
    'sEx',zeros(NElem(1,1),1),...
    'sEy',zeros(NElem(1,1),1),...
    'sEz',zeros(NElem(1,1),1),...
    'ECtoR',zeros(NElem(1,1),1),...
    'EAreaR',zeros(NElem(1,1),1),...
    'BIndS',0.0,...
    'EIndS',0.0,...
    'BIndE',0.0,...
    'EIndE',0.0);

cactusGeom.strut = repmat(strut,1,NStrut(1));

for i=1:NStrut(1,1)

    if i ~= 1
        myfgetl(fid); # skip line
        myfgetl(fid); # skip line
    end


    cactusGeom.strut(1,i).NElem  = NElem(1,1);
    TtoC = processLine(fid);
    cactusGeom.strut(1,i).TtoC = TtoC(1,1);

    MCx =  processLine(fid);
    cactusGeom.strut(1,i).MCx(1:NElem(1,1)+1,1) = MCx(1:NElem(1,1)+1,1);
    MCy =  processLine(fid);
    cactusGeom.strut(1,i).MCy(1:NElem(1,1)+1,1) = MCy(1:NElem(1,1)+1,1);
    MCz =  processLine(fid);
    cactusGeom.strut(1,i).MCz(1:NElem(1,1)+1,1) = MCz(1:NElem(1,1)+1,1);

    CtoR =  processLine(fid);
    cactusGeom.strut(1,i).CtoR(1:NElem(1,1)+1,1) = CtoR(1:NElem(1,1)+1,1);

    PEx =  processLine(fid);
    cactusGeom.strut(1,i).PEx(1:NElem(1,1),1) = PEx(1:NElem(1,1),1);
    PEy =  processLine(fid);
    cactusGeom.strut(1,i).PEy(1:NElem(1,1),1) = PEy(1:NElem(1,1),1);
    PEz =  processLine(fid);
    cactusGeom.strut(1,i).PEz(1:NElem(1,1),1) = PEz(1:NElem(1,1),1);

    sEx =  processLine(fid);
    cactusGeom.strut(1,i).sEx(1:NElem(1,1),1) = sEx(1:NElem(1,1),1);
    sEy =  processLine(fid);
    cactusGeom.strut(1,i).sEy(1:NElem(1,1),1) = sEy(1:NElem(1,1),1);
    sEz =  processLine(fid);
    cactusGeom.strut(1,i).sEz(1:NElem(1,1),1) = sEz(1:NElem(1,1),1);

    ECtoR =  processLine(fid);
    cactusGeom.strut(1,i).ECtoR(1:NElem(1,1),1) = ECtoR(1:NElem(1,1),1);
    EAreaR =  processLine(fid);
    cactusGeom.strut(1,i).EAreaR(1:NElem(1,1),1) = EAreaR(1:NElem(1,1),1);
    BIndS =  processLine(fid);
    cactusGeom.strut(1,i).BIndS = BIndS(1,1);
    EIndS =  processLine(fid);
    cactusGeom.strut(1,i).EIndS = EIndS(1,1);
    BIndE =  processLine(fid);
    cactusGeom.strut(1,i).BIndE = BIndE(1,1);
    EIndE =  processLine(fid);
    cactusGeom.strut(1,i).EIndE = EIndE(1,1);
end

fclose(fid);

end

function [data] = processLine(fid)

line = myfgetl(fid);

# Find where all of the delimiters are
delimiter_log_idx = line == ' ';
# Reduce delimiter index to just two consecutive spaces
delimiter_idx_size = 0; #Figure out the delimiter array size
for i = 2:length(delimiter_log_idx)-1
    if (delimiter_log_idx(i-1) == 1) && (delimiter_log_idx(i) == 1) && (delimiter_log_idx(i+1) == 0)
        delimiter_idx_size = delimiter_idx_size + 1;
    end
end

delimiter_idx = zeros(1,delimiter_idx_size); #Initialize variable scope
j = 1;
for i = 2:length(delimiter_log_idx)-1
    if (delimiter_log_idx(i-1) == 1) && (delimiter_log_idx(i) == 1) && (delimiter_log_idx(i+1) == 0)
        delimiter_idx(j) = i;
        j = j+1;
    end
end

delimiter_idx = [delimiter_idx,length(line)+1]; #They all have a prefix in the .geom file, so don't include a 0 index and thus don't grab the prefix.  The logic above also skips the type and blade lines, so no need to parse for those either.
data = zeros(length(delimiter_idx)-1,1);
# Extract the data from the beginning to the last delimiter
for k = 2:length(delimiter_idx)
    data(k-1) = real(str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1)));
end

end
