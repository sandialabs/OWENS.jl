function [ mac ] = calculateMacTable(file1,file2)
#calculates a modal assurance criterion table from eigenvectors


vec1 = load(file1);
vec2 = load(file2);
vec1 = vec1.eigVec;
vec2 = vec2.eigVec;

vec1 = vec1(:,2:2:end);
vec2 = vec2(:,2:2:end);

[numDof,numModes] = size(vec1);
mac = zeros(numModes);

for i=1:numModes
    for j=1:numModes
        a = vec1(:,i);
        b = vec2(:,j);
        
        mac(i,j) =  ((a'*b)*(b'*a))/((a'*a)*(b'*b));
    end
end


end
