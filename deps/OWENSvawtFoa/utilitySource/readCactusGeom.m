function [cactusGeom] = readCactusGeom(geomfn)

fid = fopen(geomfn,'r');

[cactusGeom.NBlade] = processLine(fid);
[cactusGeom.NStrut] = processLine(fid);
[cactusGeom.RotN] = processLine(fid);
[cactusGeom.RotP] = processLine(fid);
[cactusGeom.RefAR] = processLine(fid);
[cactusGeom.RefR] = processLine(fid);

cactusGeom.blade = struct('NElem',cell(1,cactusGeom.NBlade),...
'FlipN',cell(1,cactusGeom.NBlade),...
'QCx',cell(1,cactusGeom.NBlade),...
'QCy',cell(1,cactusGeom.NBlade),...
'QCz',cell(1,cactusGeom.NBlade),...
'tx',cell(1,cactusGeom.NBlade),...
'ty',cell(1,cactusGeom.NBlade),...
'tz',cell(1,cactusGeom.NBlade),...
'CtoR',cell(1,cactusGeom.NBlade),...
'PEx',cell(1,cactusGeom.NBlade),...
'PEy',cell(1,cactusGeom.NBlade),...
'PEz',cell(1,cactusGeom.NBlade),...
'tEx',cell(1,cactusGeom.NBlade),...
'tEy',cell(1,cactusGeom.NBlade),...
'tEz',cell(1,cactusGeom.NBlade),...
'nEx',cell(1,cactusGeom.NBlade),...
'nEy',cell(1,cactusGeom.NBlade),...
'nEz',cell(1,cactusGeom.NBlade),...
'sEx',cell(1,cactusGeom.NBlade),...
'sEy',cell(1,cactusGeom.NBlade),...
'sEz',cell(1,cactusGeom.NBlade),...
'ECtoR',cell(1,cactusGeom.NBlade),...
'EAreaR',cell(1,cactusGeom.NBlade),...
'iSect',cell(1,cactusGeom.NBlade));

cactusGeom.strut = struct('NElem',cell(1,cactusGeom.NStrut),...
'TtoC',cell(1,cactusGeom.NStrut),...
'MCx',cell(1,cactusGeom.NStrut),...
'MCy',cell(1,cactusGeom.NStrut),...
'MCz',cell(1,cactusGeom.NStrut),...
'CtoR',cell(1,cactusGeom.NStrut),...
'PEx',cell(1,cactusGeom.NStrut),...
'PEy',cell(1,cactusGeom.NStrut),...
'PEz',cell(1,cactusGeom.NStrut),...
'sEx',cell(1,cactusGeom.NStrut),...
'sEy',cell(1,cactusGeom.NStrut),...
'sEz',cell(1,cactusGeom.NStrut),...
'ECtoR',cell(1,cactusGeom.NStrut),...
'EAreaR',cell(1,cactusGeom.NStrut),...
'BIndS',cell(1,cactusGeom.NStrut),...
'EIndS',cell(1,cactusGeom.NStrut),...
'BIndE',cell(1,cactusGeom.NStrut),...
'EIndE',cell(1,cactusGeom.NStrut));

myfgetl(fid);

for i=1:cactusGeom.NBlade
    [cactusGeom.blade(i)] = readBladeBlock(fid);
end

for i=1:cactusGeom.NStrut
    [cactusGeom.strut(i)] = readStrutBlock(fid);
end

fclose(fid);

end

function [blade] = readBladeBlock(fid)
myfgetl(fid);
[blade.NElem] =  processLine(fid);
[blade.FlipN] =  processLine(fid);

%     fsecstr = [];
%     isecstr = [];
%     for j=1:blade.NElem+1
%         fsecstr = [fsecstr,'%f'];
%         isecstr = [isecstr,'%i'];
%     end
%     fsecel = fsecstr(1:end-2);
%     isecel = isecstr(1:end-2);

[blade.QCx] =  processLine(fid);
[blade.QCy] =  processLine(fid);
[blade.QCz] =  processLine(fid);

[blade.tx] =  processLine(fid);
[blade.ty] =  processLine(fid);
[blade.tz] =  processLine(fid);

[blade.CtoR] =  processLine(fid);

[blade.PEx] =  processLine(fid);
[blade.PEy] =  processLine(fid);
[blade.PEz] =  processLine(fid);

[blade.tEx] =  processLine(fid);
[blade.tEy] =  processLine(fid);
[blade.tEz] =  processLine(fid);

[blade.nEx] =  processLine(fid);
[blade.nEy] =  processLine(fid);
[blade.nEz] =  processLine(fid);

[blade.sEx] =  processLine(fid);
[blade.sEy] =  processLine(fid);
[blade.sEz] =  processLine(fid);

[blade.ECtoR] =  processLine(fid);
[blade.EAreaR] =  processLine(fid);
[blade.iSect] =  processLine(fid);
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
delimiter_idx = length(delimiter_log_idx); %Initialize variable scope if there are no delimiters
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
