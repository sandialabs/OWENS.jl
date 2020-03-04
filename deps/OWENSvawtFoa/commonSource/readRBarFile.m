function [joint] = readRBarFile(filename,joint,mesh)
    fid = fopen(filename,'r');
    if(fid ~= -1)
        numRbar = fscanf(fid,'%i',1);
        masterNode = zeros(numRbar);
        slaveNode = zeros(numRbar);
        lengthVec = zeros(numRbar,3);
        if(isempty(joint))
            indexStart = 0;
            joint = zeros(numRbar,8);
        else
            indexStart = joint(end,1);
        end
        for i=1:numRbar
           masterNode(i) = fscanf(fid,'%i',1);
           slaveNode(i) = fscanf(fid,'%i',1);
           lengthVec(i,1) = mesh.x(slaveNode(i)) - mesh.x(masterNode(i));
           lengthVec(i,2) = mesh.y(slaveNode(i)) - mesh.y(masterNode(i));
           lengthVec(i,3) = mesh.z(slaveNode(i)) - mesh.z(masterNode(i));

           joint(indexStart + i,1) = indexStart + i;
           joint(indexStart + i,2) = masterNode(i);
           joint(indexStart + i,3) = slaveNode(i);
           joint(indexStart + i,4) = 5; %this is the ID for a rbar constraint
           joint(indexStart + i,5:7) = lengthVec(i,:);
        end
    end

end