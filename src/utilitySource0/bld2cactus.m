function bld2cactus(bldfn,outname,normalizedBladeRef,plotFlag)

a = importdata(bldfn);

bladeNum = a(:,1);

numBlades = max(bladeNum);
numStruts = min(bladeNum);
if(numStruts>0)
    numStruts = 0;
else 
    numStruts = abs(numStruts);
end

strutStartIndex = 0;
for i=1:length(bladeNum)
    if(isnan(a(i,end)))
        strutStartIndex = i;
        break;
    end
end



if(strutStartIndex!=0)
    strutDataBlock = a(strutStartIndex:end,:);
    [strutEntries, ~] = size(strutDataBlock);
    numNodesPerStrut = strutEntries/numStruts;
    numElPerStrut = numNodesPerStrut - 1;
else
    [temp,~]=size(a);
    strutStartIndex = temp + 1;
end

bladeDataBlock = a(1:strutStartIndex-1,:);
[bladeEntries, ~] = size(bladeDataBlock);
numNodesPerBlade = bladeEntries/numBlades;
numElPerBlade = numNodesPerBlade - 1;


#calculate maximum radius of turbine
maxRadius = 0.0;
index = 1;
for i=1:numBlades
    for j=1:numNodesPerBlade
        radius_ij = sqrt(bladeDataBlock(index,5)^2 + bladeDataBlock(index,6)^2);
        if(radius_ij > maxRadius)
            maxRadius = radius_ij;
        end
        index = index + 1;
    end
end
refR = maxRadius;

nodeNum = bladeDataBlock(:,3);
elementNum = bladeDataBlock(:,4);
qc= bladeDataBlock(:,5:7)/refR; #normalized by reference radius
nvec = bladeDataBlock(:,8:10);
tvec = bladeDataBlock(:,11:13);
chord = bladeDataBlock(:,14)/refR; #normalized by reference radius
airfoil = bladeDataBlock(:,15);

index = 1;

for i=1:numBlades
    
    for j=1:numNodesPerBlade
        
        blade{i}.QCx(j) = qc(index,2); #cactus x is vawtgen y
        blade{i}.QCy(j) = qc(index,3); #cactus y is vawtgen z
        blade{i}.QCz(j) = qc(index,1); #cactus z is vawtgen x 
        
        blade{i}.nx(j) = nvec(index,2); #cactus x is vawtgen y
        blade{i}.ny(j) = nvec(index,3); #cactus y is vawtgen z
        blade{i}.nz(j) = nvec(index,1); #cactus z is vawtgen x 
        
        blade{i}.tx(j) = tvec(index,2); #cactus x is vawtgen y
        blade{i}.ty(j) = tvec(index,3); #cactus y is vawtgen z
        blade{i}.tz(j) = tvec(index,1); #cactus z is vawtgen x 
        
        blade{i}.c(j) = chord(index);
        
        if(j<numNodesPerBlade)
            elDelta = [qc(index+1,1);qc(index+1,2);qc(index+1,3)] - [qc(index,1);qc(index,2);qc(index,3)];
            elLen = norm(elDelta);
            blade{i}.area(j) = 0.5*(chord(index) + chord(index+1))*elLen;
            blade{i}.afNum(j) = airfoil(index);
        end
        
        index = index + 1;
    end
    
    if(normalizedBladeRef!=-1.0)
       clear blade{i}.area;
    #set up normalized blade ref dimension for mappiing to aero refinement
    lenx=abs(blade{i}.QCx(end)-blade{i}.QCx(1));
    leny=abs(blade{i}.QCy(end)-blade{i}.QCy(1));
    lenz=abs(blade{i}.QCz(end)-blade{i}.QCz(1));
    lenarray=[lenx;leny;lenz];
    [~,maxindex] = max(lenarray);
    bladeRefLength = lenarray(maxindex);

    for j=1:numNodesPerBlade
#         bladeFrac(j) = abs((blade{i}.QCx(j)-blade{i}.QCx(1))^2 + (blade{i}.QCy(j)-blade{i}.QCy(1))^2 + (blade{i}.QCz(j)-blade{i}.QCz(1))^2)/bladeRefLength;
          if(maxindex==1)
            bladeFrac(j) = abs((blade{i}.QCx(j)-blade{i}.QCx(1)))/bladeRefLength;
          end
          if(maxindex==2)
            bladeFrac(j) = abs((blade{i}.QCy(j)-blade{i}.QCy(1)))/bladeRefLength;
          end
          if(maxindex==3)
            bladeFrac(j) = abs((blade{i}.QCz(j)-blade{i}.QCz(1)))/bladeRefLength;
          end
            
    end
        bladeFracMid = 0.5*(bladeFrac(1:end-1) + bladeFrac(2:end));

        normBladeRefMid = 0.5*(normalizedBladeRef(1:end-1) + normalizedBladeRef(2:end));
    
    blade{i}.QCx = interp1(bladeFrac,blade{i}.QCx,normalizedBladeRef,'linear');
    blade{i}.QCy = interp1(bladeFrac,blade{i}.QCy,normalizedBladeRef,'linear');
    blade{i}.QCz = interp1(bladeFrac,blade{i}.QCz,normalizedBladeRef,'linear');
    
    blade{i}.nx = interp1(bladeFrac,blade{i}.nx,normalizedBladeRef,'linear');
    blade{i}.ny = interp1(bladeFrac,blade{i}.ny,normalizedBladeRef,'linear');
    blade{i}.nz = interp1(bladeFrac,blade{i}.nz,normalizedBladeRef,'linear');
    
    blade{i}.tx = interp1(bladeFrac,blade{i}.tx,normalizedBladeRef,'linear');
    blade{i}.ty = interp1(bladeFrac,blade{i}.ty,normalizedBladeRef,'linear');
    blade{i}.tz = interp1(bladeFrac,blade{i}.tz,normalizedBladeRef,'linear');
    
    blade{i}.c = interp1(bladeFrac,blade{i}.c,normalizedBladeRef,'linear');
    
    blade{i}.afNum = interp1(bladeFracMid,blade{i}.afNum,normBladeRefMid,'nearest');
    
    for j=1:length(normalizedBladeRef)-1
            elDelta = [blade{i}.QCx(j+1);blade{i}.QCy(j+1);blade{i}.QCz(j+1)] - [blade{i}.QCx(j);blade{i}.QCy(j);blade{i}.QCz(j)];
            elLen = norm(elDelta);
            blade{i}.area(j) = 0.5*(blade{i}.c(j+1) + blade{i}.c(j))*elLen;
    end
        blade{i}.area = blade{i}.area/(refR^2); #normalize by refAr
    end
end


index = 0;
strut=[];
for i=1:numStruts
    for j=1:numNodesPerStrut
        index = index + 1;
        strut{i}.SEx(j) = strutDataBlock(index,5)/refR;
        strut{i}.SEy(j) = strutDataBlock(index,6)/refR;
        strut{i}.SEz(j) = strutDataBlock(index,7)/refR;
        strut{i}.c(j) = strutDataBlock(index,8)/refR;
        strut{i}.tc(j) = strutDataBlock(index,9);
        strut{i}.bladeIndex(j) = strutDataBlock(index,10);
        
        if(j<numNodesPerStrut)
            elDelta = [strutDataBlock(index+1,5);strutDataBlock(index+1,6);strutDataBlock(index+1,7)] - [strutDataBlock(index,5);strutDataBlock(index,6);strutDataBlock(index,7)];
            elLen = norm(elDelta);
            strut{i}.area(j) = 0.5*(strutDataBlock(index,8) + strutDataBlock(index+1,8))*elLen;
        end
        
    end
        strut{i}.area = strut{i}.area/(refR^2); #normalize by refAr
end

#search for coincident strut node and blade node
for i=1:numStruts
    strutEndCoord = [strut{i}.SEx(end); strut{i}.SEy(end); strut{i}.SEz(end)];
    bindex=strut{i}.bladeIndex(1);
    for j=1:numNodesPerBlade
        delta = [blade{bindex}.QCx(j); blade{bindex}.QCy(j); blade{bindex}.QCz(j)]-strutEndCoord;
        if(norm(delta)<1.0e-4)
            bladeNodeIndex = j;
            strut{i}.bladeElConnectionIndex = bladeNodeIndex;
            break;
        end
    end
end


#calculate projected frontal area of volume swept by turbine (assumed all
#blades are the same
A = 0;
for i=1:(length(blade{1}.QCx)-1)
    
    dely = blade{1}.QCy(i+1) - blade{1}.QCy(i);
    delA = 0.5*(blade{1}.QCz(i+1) + blade{1}.QCz(i))*dely;
    A=A+abs(delA);
    
end
refAr = 2.0*A;

#may want to correct for quarter chord vs. semi chord of the strut here, 
#before writing

writeCactus(outname,blade,strut,refAr,refR);

if(plotFlag)
    orientationVectorScaleFactor = 0.05*max(max(abs(qc)));
    plotCactusGeom(blade,strut,orientationVectorScaleFactor);
end



end

function writeCactus(outname,blade,strut,refAr,refR)
    
    fid=fopen(outname,'w');
    writeCactusHeader(fid,length(blade),length(strut),refAr,refR);
    
    for i=1:length(blade)
        writeCactusBlade(fid,blade{i},i)
    end
    
    for i=1:length(strut)
        writeCactusStrut(fid,strut{i},i)
    end
    
    fclose(fid);

end

function writeCactusHeader(fid,numBlades,numStruts,refAr,refR)
    fprintf(fid,'NBlade: #3d \n',numBlades);
    fprintf(fid,'NStrut: #3d \n',numStruts);
    fprintf(fid,'RotN: '); printScientific(fid,[0,1,0]); fprintf(fid,'\n');
    fprintf(fid,'RotP: '); printScientific(fid,[0,0,0]); fprintf(fid,'\n');
    fprintf(fid,'RefAR: '); printScientific(fid,refAr); fprintf(fid,'\n');
    fprintf(fid,'RefR: '); printScientific(fid,refR); fprintf(fid,'\n');
    fprintf(fid,'Type: VAWT\n');
end

function writeCactusBlade(fid,blade,i)
    fprintf(fid,['Blade ',num2str(i),': \n']);
    fprintf(fid,'\tNElem: #3d \n',length(blade.afNum));
#     fprintf(fid,
    fprintf(fid,'\tQCx: ');printScientific(fid,blade.QCx); fprintf(fid,'\n');
    fprintf(fid,'\tQCy: ');printScientific(fid,blade.QCy); fprintf(fid,'\n');
    fprintf(fid,'\tQCz: ');printScientific(fid,blade.QCz); fprintf(fid,'\n');
    fprintf(fid,'\tnx: ');printScientific(fid,blade.nx); fprintf(fid,'\n');
    fprintf(fid,'\tny: ');printScientific(fid,blade.ny); fprintf(fid,'\n');
    fprintf(fid,'\tnz: ');printScientific(fid,blade.nz); fprintf(fid,'\n');
    fprintf(fid,'\ttx: ');printScientific(fid,blade.tx); fprintf(fid,'\n');
    fprintf(fid,'\tty: ');printScientific(fid,blade.ty); fprintf(fid,'\n');
    fprintf(fid,'\ttz: ');printScientific(fid,blade.tz); fprintf(fid,'\n');
    fprintf(fid,'\tCtoR: ');printScientific(fid,blade.c); fprintf(fid,'\n');
    fprintf(fid,'\tAreaR: ');printScientific(fid,blade.area); fprintf(fid,'\n');
    fprintf(fid,'\tiSect: ');fprintf(fid,'#3d ',blade.afNum); fprintf(fid,'\n');
end

function writeCactusStrut(fid,strut,i)
    fprintf(fid,['Strut ',num2str(i),': \n']);
    fprintf(fid,'\tNElem: #3d \n',length(strut.area));
    fprintf(fid,'\tSEx: ');printScientific(fid,strut.SEx); fprintf(fid,'\n');
    fprintf(fid,'\tSEy: ');printScientific(fid,strut.SEy); fprintf(fid,'\n');
    fprintf(fid,'\tSEz: ');printScientific(fid,strut.SEz); fprintf(fid,'\n');
    fprintf(fid,'\tCtoR: ');printScientific(fid,strut.c); fprintf(fid,'\n');
    fprintf(fid,'\tAreaR: ');printScientific(fid,strut.area); fprintf(fid,'\n');
    fprintf(fid,'\tTtoC: ');printScientific(fid,strut.tc(1)); fprintf(fid,'\n');
    fprintf(fid,'\tBInd: #3d \n',strut.bladeIndex(1));
    fprintf(fid,'\tIEnd: #3d \n',strut.bladeElConnectionIndex);
end

function plotCactusGeom(blade,strut,sf)
	close all;
    for i=1:length(blade)
        plotBlade(blade{i},sf);
    end
    for i=1:length(strut)
        plotStrut(strut{i})
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

function plotBlade(blade,sf)
    hold on;
    plot3(blade.QCx,blade.QCy,blade.QCz,'r-o');
    view(3);
    axis equal;
    
    #plot norm vecs
    for j=1:length(blade.QCx)
        a0=[blade.QCx(j),blade.QCy(j),blade.QCz(j)];
        a1=[blade.QCx(j),blade.QCy(j),blade.QCz(j)] + [blade.nx(j),blade.ny(j),blade.nz(j)].*sf;
        b1=[blade.QCx(j),blade.QCy(j),blade.QCz(j)] + [blade.tx(j),blade.ty(j),blade.tz(j)].*sf;
        plot3([a0(1) a1(1)],[a0(2) a1(2)],[a0(3) a1(3)],'-k');
        plot3([a0(1) b1(1)],[a0(2) b1(2)],[a0(3) b1(3)],'-g');
    end
end

function plotStrut(strut)
    hold on;
    plot3(strut.SEx,strut.SEy,strut.SEz,'b-o');
    view(3);
    axis equal;
end

function printScientific(fid,a)
    A_str = sprintf('#6.5e ',a);
    A_str = strrep(A_str, 'e+0','e+');
    A_str = strrep(A_str, 'e-0','e-');
    
    fprintf(fid,A_str);
end
