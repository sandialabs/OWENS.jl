function [phase1,phase2] = readResults(resultsFile,modeNum,numNodes)

fid = fopen(resultsFile,'r');

modecount = 1;

while(modecount <= modeNum && ~feof(fid))
    dum = fscanf(fid,'%s',17);
    for j=1:numNodes
       temp1 = fscanf(fid,'%e',6);
       phase1(j,:,modecount) = temp1';
    end

    dum = fscanf(fid,'%s',10);
    for j=1:numNodes
       temp2 = fscanf(fid,'%e',6);
       phase2(j,:,modecount) = temp2';
    end
    modecount = modecount + 1;
end
fclose(fid);

if(modeNum > modecount)
   error('modeNum is greater than number of modes in results file. Exiting');
end

phase1 = phase1(:,:,modeNum);
phase2 = phase2(:,:,modeNum);

end

