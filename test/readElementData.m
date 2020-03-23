function [el] = readElementData(numElements,elfile,ortfile,bladeData_struct)
%readElementData  reads element data
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [el] = readElementData(numElements,elfile,ortfile,bldfile
%
%   This function reads element data and stores data in the element data
%   object.
%
%      input:
%      numElements   = number of elements in structural mesh
%      elfile        = element data filename
%      ortfile       = element orientation filename
%      bldfile       = blade data filename

%      output:
%      el            = element data object

fid = fopen(elfile,'r'); %open element data file

temp_sectionPropsArray = struct('ac',zeros(1,2),...
'twist',zeros(1,2),...
'rhoA',zeros(1,2),...
'EIyy',zeros(1,2),...
'EIzz',zeros(1,2),...
'GJ',zeros(1,2),...
'EA',zeros(1,2),...
'rhoIyy',zeros(1,2),...
'rhoIzz',zeros(1,2),...
'rhoJ',zeros(1,2),...
'zcm',zeros(1,2),...
'ycm',zeros(1,2),...
'a',zeros(1,2),...
'EIyz',zeros(1,2),...
'alpha1',zeros(1,2),...
'alpha2',zeros(1,2),...
'alpha3',zeros(1,2),...
'alpha4',zeros(1,2),...
'alpha5',zeros(1,2),...
'alpha6',zeros(1,2),...
'rhoIyz',zeros(1,2),...
'b',zeros(1,2),...
'a0',zeros(1,2),...
'aeroCenterOffset',zeros(1,2));

sectionPropsArray = repmat(temp_sectionPropsArray,1,numElements);

data1 = zeros(1,17);
data2 = zeros(1,17);
for i=1:numElements
    data1=getSplitLine(fid,'	'); %read element data
    data2=getSplitLine(fid,'	');

    %structural properties
    sectionPropsArray(i).ac = -([data1(2) data2(2)]-0.5);
    sectionPropsArray(i).twist=[data1(3) data2(3)];
    sectionPropsArray(i).rhoA = [data1(4) data2(4)];
    sectionPropsArray(i).EIyy = [data1(5) data2(5)];
    sectionPropsArray(i).EIzz = [data1(6) data2(6)];
    if(abs(sectionPropsArray(i).EIyy - sectionPropsArray(i).EIzz) < 1.0e-3)
        sectionPropsArray(i).EIzz = sectionPropsArray(i).EIzz*1.0001;
    end
    sectionPropsArray(i).GJ = [data1(7) data2(7)];
    sectionPropsArray(i).EA = [data1(8) data2(8)];
    sectionPropsArray(i).alpha1 = [data1(9) data2(9)];

    sectionPropsArray(i).rhoIyy = [data1(10) data2(10)];
    sectionPropsArray(i).rhoIzz = [data1(11) data2(11)];
    sectionPropsArray(i).rhoJ = [(data1(10)+data1(11)), (data2(10)+data2(11))];
    sectionPropsArray(i).zcm = [data1(14) data2(14)];
    sectionPropsArray(i).ycm = [data1(15) data2(15)];
    sectionPropsArray(i).a = [data1(17) data2(17)];

    %coupling factors
    sectionPropsArray(i).EIyz = [0 0];
    sectionPropsArray(i).alpha1 = [0 0];
    sectionPropsArray(i).alpha2 = [0 0];
    sectionPropsArray(i).alpha3 = [0 0];
    sectionPropsArray(i).alpha4 = [0 0];
    sectionPropsArray(i).alpha5 = [0 0];
    sectionPropsArray(i).alpha6 = [0 0];
    sectionPropsArray(i).rhoIyz = [0 0];
    sectionPropsArray(i).b = [0 0];
    sectionPropsArray(i).a0 = [2*pi 2*pi];

end

nodeNum = bladeData_struct.nodeNum;  %node number associated with blade section
elNum = bladeData_struct.elementNum;    %element number associated with blade sectino
bladeData = bladeData_struct.remaining;  %blade data

chord = zeros(length(elNum),1);
for i=1:length(elNum)
    chord(nodeNum(i)) = bladeData(i,10);  %store chord of blade sections
end

for i=1:length(elNum)
    if (elNum(i)~=-1)

        sectionPropsArray(elNum(i)).b = 0.5.*[chord(nodeNum(i)) chord(nodeNum(i+1))]; %element semi chord
        sectionPropsArray(elNum(i)).a0 = [bladeData(i,12) bladeData(i+1,12)];         %element lift curve slope (needed for flutter analysis)

        %convert "a" to semichord fraction aft of halfchord
        sectionPropsArray(elNum(i)).a = (sectionPropsArray(elNum(i)).a + 0.25*(2*sectionPropsArray(elNum(i)).b) - sectionPropsArray(elNum(i)).b)./sectionPropsArray(elNum(i)).b;

        %convert "ac" to semichord fraction foreward of halfchord
        sectionPropsArray(elNum(i)).ac = (sectionPropsArray(elNum(i)).ac).*2;

        %physical aero center offset from elastic axis
        sectionPropsArray(elNum(i)).aeroCenterOffset = (sectionPropsArray(elNum(i)).ac).*sectionPropsArray(elNum(i)).b - sectionPropsArray(elNum(i)).a;
    end
end


disp('EIyz, rhoIyz deactivated');
fclose(fid); %close element file

%read element orientation data
elLen = zeros(numElements,1);
psi = zeros(numElements,1);
theta = zeros(numElements,1);
roll = zeros(numElements,1);
fid = fopen(ortfile,'r');
temp = zeros(1,5);
for i=1:numElements
    temp = getSplitLine(fid,'	');
    elLen(i)=temp(5);
    psi(i)=temp(2);
    theta(i)=temp(3);
    roll(i)=temp(4);
end

%store data in element object
el.props = sectionPropsArray;
el.elLen=elLen;
el.psi=psi;
el.theta=theta;
el.roll=roll;

el.rotationalEffects = true(1,numElements);

fclose(fid); %close ort file
end
