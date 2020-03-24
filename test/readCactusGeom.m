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

blade = struct('NElem',0.0,...
    'FlipN',0.0,...
    'QCx',0.0,...
    'QCy',0.0,...
    'QCz',0.0,...
    'tx',0.0,...
    'ty',0.0,...
    'tz',0.0,...
    'CtoR',0.0,...
    'PEx',0.0,...
    'PEy',0.0,...
    'PEz',0.0,...
    'tEx',0.0,...
    'tEy',0.0,...
    'tEz',0.0,...
    'nEx',0.0,...
    'nEy',0.0,...
    'nEz',0.0,...
    'sEx',0.0,...
    'sEy',0.0,...
    'sEz',0.0,...
    'ECtoR',0.0,...
    'EAreaR',0.0,...
    'iSect',0.0);

cactusGeom.blade = repmat(blade,1,NBlade(1));

strut = struct('NElem',0.0,...
    'TtoC',0.0,...
    'MCx',0.0,...
    'MCy',0.0,...
    'MCz',0.0,...
    'CtoR',0.0,...
    'PEx',0.0,...
    'PEy',0.0,...
    'PEz',0.0,...
    'sEx',0.0,...
    'sEy',0.0,...
    'sEz',0.0,...
    'ECtoR',0.0,...
    'EAreaR',0.0,...
    'BIndS',0.0,...
    'EIndS',0.0,...
    'BIndE',0.0,...
    'EIndE',0.0);

cactusGeom.strut = repmat(strut,1,NStrut(1));

myfgetl(fid);

for i=1:cactusGeom.NBlade

    myfgetl(fid);
    cactusGeom.blade(1,i).NElem =  processLine(fid);
    cactusGeom.blade(1,i).FlipN =  processLine(fid);

    cactusGeom.blade(1,i).QCx =  processLine(fid);
    cactusGeom.blade(1,i).QCy =  processLine(fid);
    cactusGeom.blade(1,i).QCz =  processLine(fid);

    cactusGeom.blade(1,i).tx =  processLine(fid);
    cactusGeom.blade(1,i).ty =  processLine(fid);
    cactusGeom.blade(1,i).tz =  processLine(fid);

    cactusGeom.blade(1,i).CtoR =  processLine(fid);

    cactusGeom.blade(1,i).PEx =  processLine(fid);
    cactusGeom.blade(1,i).PEy =  processLine(fid);
    cactusGeom.blade(1,i).PEz =  processLine(fid);

    cactusGeom.blade(1,i).tEx =  processLine(fid);
    cactusGeom.blade(1,i).tEy =  processLine(fid);
    cactusGeom.blade(1,i).tEz =  processLine(fid);

    cactusGeom.blade(1,i).nEx =  processLine(fid);
    cactusGeom.blade(1,i).nEy =  processLine(fid);
    cactusGeom.blade(1,i).nEz =  processLine(fid);

    cactusGeom.blade(1,i).sEx =  processLine(fid);
    cactusGeom.blade(1,i).sEy =  processLine(fid);
    cactusGeom.blade(1,i).sEz =  processLine(fid);

    cactusGeom.blade(1,i).ECtoR =  processLine(fid);
    cactusGeom.blade(1,i).EAreaR =  processLine(fid);
    cactusGeom.blade(1,i).iSect =  processLine(fid);
end

for i=1:cactusGeom.NStrut
    [Strut] = readStrutBlock(fid);
    [cactusGeom.strut(i)]  = [Strut];
end

fclose(fid);

end

function [blade] = readStrutBlock(fid)
myfgetl(fid);
[blade.NElem] =  processLine(fid);

%     fsecstr = [];
%     isecstr = [];
%     for j=1:blade.NElem+1
%         fsecstr = [fsecstr,'%f'];
%         isecstr = [isecstr,'%i'];
%     end
%     fsecel = fsecstr(1:end-2);
%     isecel = isecstr(1:end-2);

[blade.TtoC] = processLine(fid);

[blade.MCx] =  processLine(fid);
[blade.MCy] =  processLine(fid);
[blade.MCz] =  processLine(fid);

[blade.CtoR] =  processLine(fid);

[blade.PEx] =  processLine(fid);
[blade.PEy] =  processLine(fid);
[blade.PEz] =  processLine(fid);

[blade.sEx] =  processLine(fid);
[blade.sEy] =  processLine(fid);
[blade.sEz] =  processLine(fid);

[blade.ECtoR] =  processLine(fid);
[blade.EAreaR] =  processLine(fid);
[blade.BIndS] =  processLine(fid);
[blade.EIndS] =  processLine(fid);
[blade.BIndE] =  processLine(fid);
[blade.EIndE] =  processLine(fid);
end

function [data] = processLine(fid)

line = myfgetl(fid);

% Find where all of the delimiters are
delimiter_log_idx = line == ' ';
% Reduce delimiter index to just two consecutive spaces
delimiter_idx_size = 0; %Figure out the delimiter array size
for i = 2:length(delimiter_log_idx)-1
    if (delimiter_log_idx(i-1) == 1) && (delimiter_log_idx(i) == 1) && (delimiter_log_idx(i+1) == 0)
        delimiter_idx_size = delimiter_idx_size + 1;
    end
end

delimiter_idx = zeros(1,delimiter_idx_size); %Initialize variable scope
j = 1;
for i = 2:length(delimiter_log_idx)-1
    if (delimiter_log_idx(i-1) == 1) && (delimiter_log_idx(i) == 1) && (delimiter_log_idx(i+1) == 0)
        delimiter_idx(j) = i;
        j = j+1;
    end
end

delimiter_idx = [delimiter_idx,length(line)+1]; %They all have a prefix in the .geom file, so don't include a 0 index and thus don't grab the prefix.  The logic above also skips the type and blade lines, so no need to parse for those either.
data = zeros(length(delimiter_idx)-1,1);
% Extract the data from the beginning to the last delimiter
for k = 2:length(delimiter_idx)
    data(k-1) = real(str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1)));
end

end
