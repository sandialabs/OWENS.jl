function [cactusGeom] = readCactusGeom(geomfn)

fid = fopen(geomfn,'r');

[cactusGeom.NBlade,fid] = processLine(fid,'','%[NBlade:]','%i');
[cactusGeom.NStrut,fid] = processLine(fid,'','%[NStrut:]','%i');
[cactusGeom.RotN,fid] = processLine(fid,'','%[RotN:]','%f%f%f');
[cactusGeom.RotP,fid] = processLine(fid,'','%[RotP:]','%f%f%f');
[cactusGeom.RefAR,fid] = processLine(fid,'','%[RefAR:]','%f');
[cactusGeom.RefR,fid] = processLine(fid,'','%[RefR:]','%f');
fgetl(fid);

for i=1:cactusGeom.NBlade
    [cactusGeom.blade(i),fid] = readBladeBlock(fid);
end

for i=1:cactusGeom.NStrut
   [cactusGeom.strut(i),fid] = readStrutBlock(fid);
end

fclose(fid);

end

function [blade,fid] = readBladeBlock(fid)
    dum=fgetl(fid);
    [blade.NElem,fid] =  processLine(fid,'\t','%[NElem:]','%i');
    [blade.FlipN,fid] =  processLine(fid,'\t','%[FlipN:]','%i');
    
    fsecstr = [];
    isecstr = [];
    for j=1:blade.NElem+1;
        fsecstr = [fsecstr,'%f'];
        isecstr = [isecstr,'%i'];
    end
        fsecel = fsecstr(1:end-2);
        isecel = isecstr(1:end-2);
    
    [blade.QCx,fid] =  processLine(fid,'\t','%[QCx:]',fsecstr);
    [blade.QCy,fid] =  processLine(fid,'\t','%[QCy:]',fsecstr);
    [blade.QCz,fid] =  processLine(fid,'\t','%[QCz:]',fsecstr);
    
    [blade.tx,fid] =  processLine(fid,'\t','%[tx:]',fsecstr);
    [blade.ty,fid] =  processLine(fid,'\t','%[ty:]',fsecstr);
    [blade.tz,fid] =  processLine(fid,'\t','%[tz:]',fsecstr);
    
    [blade.CtoR,fid] =  processLine(fid,'\t','%[CtoR:]',fsecstr);
    
    [blade.PEx,fid] =  processLine(fid,'\t','%[PEx:]',fsecel);
    [blade.PEy,fid] =  processLine(fid,'\t','%[PEy:]',fsecel);
    [blade.PEz,fid] =  processLine(fid,'\t','%[PEz:]',fsecel);
    
    [blade.tEx,fid] =  processLine(fid,'\t','%[tEx:]',fsecel);
    [blade.tEy,fid] =  processLine(fid,'\t','%[tEy:]',fsecel);
    [blade.tEz,fid] =  processLine(fid,'\t','%[tEz:]',fsecel);
    
    [blade.nEx,fid] =  processLine(fid,'\t','%[nEx:]',fsecel);
    [blade.nEy,fid] =  processLine(fid,'\t','%[nEy:]',fsecel);
    [blade.nEz,fid] =  processLine(fid,'\t','%[nEz:]',fsecel);
    
    [blade.sEx,fid] =  processLine(fid,'\t','%[sEx:]',fsecel);
    [blade.sEy,fid] =  processLine(fid,'\t','%[sEy:]',fsecel);
    [blade.sEz,fid] =  processLine(fid,'\t','%[sEz:]',fsecel);
    
    [blade.ECtoR,fid] =  processLine(fid,'\t','%[ECtoR:]',fsecel);
    [blade.EAreaR,fid] =  processLine(fid,'\t','%[EAreaR:]',fsecel);
    [blade.iSect,fid] =  processLine(fid,'\t','%[iSect:]',isecel);
end

function [blade,fid] = readStrutBlock(fid)
    dum=fgetl(fid);
    [blade.NElem,fid] =  processLine(fid,'\t','%[NElem:]','%i');
        
    fsecstr = [];
    isecstr = [];
    for j=1:blade.NElem+1;
        fsecstr = [fsecstr,'%f'];
        isecstr = [isecstr,'%i'];
    end
        fsecel = fsecstr(1:end-2);
        isecel = isecstr(1:end-2);
    
    [blade.TtoC,fid] = processLine(fid,'\t','%[TtoC:]','%f');
    
    [blade.MCx,fid] =  processLine(fid,'\t','%[MCx:]',fsecstr);
    [blade.MCy,fid] =  processLine(fid,'\t','%[MCy:]',fsecstr);
    [blade.MCz,fid] =  processLine(fid,'\t','%[MCz:]',fsecstr);
    
    [blade.CtoR,fid] =  processLine(fid,'\t','%[CtoR:]',fsecstr);
    
    [blade.PEx,fid] =  processLine(fid,'\t','%[PEx:]',fsecel);
    [blade.PEy,fid] =  processLine(fid,'\t','%[PEy:]',fsecel);
    [blade.PEz,fid] =  processLine(fid,'\t','%[PEz:]',fsecel);
    
    [blade.sEx,fid] =  processLine(fid,'\t','%[sEx:]',fsecel);
    [blade.sEy,fid] =  processLine(fid,'\t','%[sEy:]',fsecel);
    [blade.sEz,fid] =  processLine(fid,'\t','%[sEz:]',fsecel);
    
    [blade.ECtoR,fid] =  processLine(fid,'\t','%[ECtoR:]',fsecel);
    [blade.EAreaR,fid] =  processLine(fid,'\t','%[EAreaR:]',fsecel);
    [blade.BIndS,fid] =  processLine(fid,'\t','%[BIndS:]','%i');
    [blade.EIndS,fid] =  processLine(fid,'\t','%[EIndS:]','%i');
    [blade.BIndE,fid] =  processLine(fid,'\t','%[BIndE:]','%i');
    [blade.EIndE,fid] =  processLine(fid,'\t','%[EIndE:]','%i');
end

function [data,fid] = processLine(fid,preFix,labelStr,numStr)
    temp = fgetl(fid);
    [data] = sscanf(temp,[preFix,labelStr,numStr]);
    len = length(labelStr)-3;

    data = data(len+1:end);

end