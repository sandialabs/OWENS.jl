function [time,ForceValHist,ForceDof] = mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)
% function [time,aeroDistLoadsArrayTime,aeroDistLoadsNodeMap,aeroDistLoadsElMap] = mapCactusLoadsFile(geomFn,loadsFn,bldFn,elFn,ortFn,meshFn)

[cactusGeom] = readCactusGeom(geomFn);
blade = cactusGeom.blade;

data = importCactusFile(loadsFn,1,2002,22,',');

%define these from params file
ft2m = 1 / 3.281;
rho = 1.225;
%     RefAR = cactusGeom.RefAR*ft2m*ft2m;
RefR = cactusGeom.RefR*ft2m;
V = 25; %m/s

normTime = data(:,1);

numAeroEl = 0;
for i=1:cactusGeom.NBlade
    numAeroEl = numAeroEl + cactusGeom.blade(i).NElem;
end

[len,~] = size(data);

numAeroTS = len/numAeroEl;
time = normTime(1:numAeroEl:end).*RefR/V;

urel = data(:,12);
uloc = urel.*V;

cn = data(:,17);
ct = data(:,18);
cm25 = data(:,15);

%     cl = data(:,13);
%     cd = data(:,14);
%
%     cx = data(:,19);
%     cy = data(:,20);
%     cz = data(:,21);

%calculate element areas
%     Fx = zeros(len,1);
%     Fy = Fx;
%     Fz = Fx;

NperSpan = zeros(len,1);
TperSpan = zeros(len,1);
M25perSpan = zeros(len,1);
Mecc = zeros(len,1);

for i=1:len

    NperSpan(i) =  cn(i)  * 0.5*rho*uloc(i)^2*(blade(data(i,3)).ECtoR(data(i,4))*RefR);
    TperSpan(i) =  ct(i)  * 0.5*rho*uloc(i)^2*(blade(data(i,3)).ECtoR(data(i,4))*RefR);
    M25perSpan(i) = cm25(i) * 0.5*rho*uloc(i)^2*(blade(data(i,3)).ECtoR(data(i,4))*RefR)*blade(data(i,3)).ECtoR(data(i,4))*RefR;
    momentArm = 0.0;
    Mecc(i) = NperSpan(i) * momentArm;
end

% Initialize bladeForce
bladeForce = struct('N',cell(1,cactusGeom.NBlade),'T',cell(1,cactusGeom.NBlade),'M25',cell(1,cactusGeom.NBlade));
for i = 1:cactusGeom.NBlade
    bladeForce(i).N = zeros(numAeroTS,blade(1).NElem);
    bladeForce(i).T = zeros(numAeroTS,blade(1).NElem);
    bladeForce(i).M25 = zeros(numAeroTS,blade(1).NElem);
end

index = 1;
for i=1:numAeroTS
    for j=1:cactusGeom.NBlade
        for k=1:blade(j).NElem
            %%
            %                 bladeForce(j).Fx(i,k) = Fx(index);
            %                 bladeForce(j).Fy(i,k) = Fy(index);
            %                 bladeForce(j).Fz(i,k) = Fz(index);
            %                 bladeForce(j).M25(i,k) = M25perSpan(index);

            %                 spanVec = [blade(j).sEx(k);blade(j).sEy(k);blade(j).sEz(k)];
            %                 tanVec = [blade(j).tEx(k);blade(j).tEy(k);blade(j).tEz(k)];
            %                 normVec = [blade(j).nEx(k);blade(j).nEy(k);blade(j).nEz(k)];
            % %                 dcm = [spanVec, tanVec, normVec];
            %                 dcm = [tanVec, spanVec, -normVec];
            %%
            bladeForce(j).N(i,k) = NperSpan(index);
            bladeForce(j).T(i,k) = TperSpan(index);
            bladeForce(j).M25(i,k) = M25perSpan(index);
            index = index + 1;
        end
    end
end

spanLocNorm = zeros(cactusGeom.NBlade,blade(1).NElem);
for i=1:cactusGeom.NBlade
    spanLocNorm(i,:) = blade(i).PEy.*RefR/(blade(i).QCy(end)*RefR);
end

[bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers] = readBldFile(bldFn);

%Initialize structuralLoad
structuralLoad = struct('N',cell(1,cactusGeom.NBlade),'T',cell(1,cactusGeom.NBlade),'M25',cell(1,cactusGeom.NBlade));
for i = 1:cactusGeom.NBlade
    structuralLoad(i).N = zeros(numAeroTS,length(structuralSpanLocNorm(1,:)));
    structuralLoad(i).T = zeros(numAeroTS,length(structuralSpanLocNorm(1,:)));
    structuralLoad(i).M25 = zeros(numAeroTS,length(structuralSpanLocNorm(1,:)));
end

for i=1:cactusGeom.NBlade
    for j=1:numAeroTS
        structuralLoad(i).N(j,:) = interp1(spanLocNorm(i,:),bladeForce(i).N(j,:),structuralSpanLocNorm(i,:),'linear','extrap');
        structuralLoad(i).T(j,:) = interp1(spanLocNorm(i,:),bladeForce(i).T(j,:),structuralSpanLocNorm(i,:),'linear','extrap');
        structuralLoad(i).M25(j,:) = interp1(spanLocNorm(i,:),bladeForce(i).M25(j,:),structuralSpanLocNorm(i,:),'linear','extrap');
    end
end

[~,numNodesPerBlade] = size(structuralNodeNumbers);

%integrate over elements

%read element data in
[mesh] = readMesh(meshFn);
[el] = readElementData(mesh.numEl,elFn,ortFn,bladeData);
numDofPerNode = 6;
%     [~,~,timeLen] = size(aeroDistLoadsArrayTime);
Fg = zeros(max(max(structuralNodeNumbers))*6,numAeroTS);
for i=1:numAeroTS
    for j = 1:cactusGeom.NBlade
        for k = 1:numNodesPerBlade-1
            %get element data
            % orientation angle,xloc,sectionProps,element order]
            elNum = structuralElNumbers(j,k);
            %get dof map
            node1 = structuralNodeNumbers(j,k);
            node2 = structuralNodeNumbers(j,k+1);
            dofList = [(node1-1)*numDofPerNode+(1:6), (node2-1)*numDofPerNode+(1:6)];

            elInput.elementOrder = 1;
            elInput.x = [mesh.x(node1), mesh.x(node2)];
            elLength = sqrt((mesh.x(node2)-mesh.x(node1))^2 + (mesh.y(node2)-mesh.y(node1))^2 + (mesh.z(node2)-mesh.z(node1))^2);
            elInput.xloc = [0 elLength];
            elInput.sectionProps.twist = el.props{elNum}.twist;
            elInput.sweepAngle = el.psi(elNum);
            elInput.coneAngle = el.theta(elNum);
            elInput.rollAngle = el.roll(elNum);

            elInput.extDistF2Node =  [structuralLoad(j).T(i,k),   structuralLoad(j).T(i,k+1)];
            elInput.extDistF3Node = -[structuralLoad(j).N(i,k),   structuralLoad(j).N(i,k+1)];
            elInput.extDistF4Node = -[structuralLoad(j).M25(i,k), structuralLoad(j).M25(i,k+1)];

            [output] = calculateLoadVecFromDistForce(elInput);
            Fe = output.Fe;

            %asssembly
            for m = 1:length(dofList)
                Fg(dofList(m),i) =  Fg(dofList(m),i)+Fe(m);
            end

        end
    end
end

%reduce Fg to nonzero components
%assumes any loaded DOF will never be identically zero throughout time
%history
ForceValHist = zeros(sum(Fg(:,1)~=0),length(Fg(1,:)));
ForceDof = zeros(sum(Fg(:,1)~=0),1);
index = 1;
for i=1:max(max(structuralNodeNumbers))*6
    if(~isempty(find(Fg(i,:), 1)))
        ForceValHist(index,:) = Fg(i,:);
        ForceDof(index) = i;
        index = index + 1;
    end
end


end

function [bladeData,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers] = readBldFile(bldFn)
%% READ IN BLD FILE
a = importCactusFile(bldFn,0,60,16,'	');

bladeNum = a(:,1);

numBlades = max(bladeNum);
% numStruts = min(bladeNum);
% if(numStruts>0)
% %     numStruts = 0;
% else
% %     numStruts = abs(numStruts);
% end

strutStartIndex = 0;
for i=1:length(bladeNum)
    if(isnan(a(i,end)))
        strutStartIndex = i;
        break;
    end
end



if(strutStartIndex~=0)
    %     strutDataBlock = a(strutStartIndex:end,:);
    %     [strutEntries, ~] = size(strutDataBlock);
    %     numNodesPerStrut = strutEntries/numStruts;
    %     numElPerStrut = numNodesPerStrut - 1;
else
    [temp,~]=size(a);
    strutStartIndex = temp + 1;
end

bladeDataBlock = a(1:strutStartIndex-1,:);
[bladeEntries, ~] = size(bladeDataBlock);
numNodesPerBlade = bladeEntries/numBlades;
% numElPerBlade = numNodesPerBlade - 1;

% [len,~]=size(bladeDataBlock);
structuralSpanLocNorm = zeros(numBlades,numNodesPerBlade);
structuralNodeNumbers = zeros(numBlades,numNodesPerBlade);
structuralElNumbers = zeros(numBlades,numNodesPerBlade);
for i=1:numBlades
    structuralSpanLocNorm(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,2)./bladeDataBlock(i*numNodesPerBlade,2);
    structuralNodeNumbers(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,3);
    structuralElNumbers(i,:) = bladeDataBlock((i-1)*numNodesPerBlade+1:1:i*numNodesPerBlade,4);
end

bladeData.numBlades = numBlades;  %assign data to bladeData object
bladeData.bladeNum = bladeDataBlock(:,1);
bladeData.h = bladeDataBlock(:,2);
bladeData.nodeNum = bladeDataBlock(:,3);
bladeData.elementNum = bladeDataBlock(:,4);
bladeData.remaining = bladeDataBlock(:,5:end);
%%
end
