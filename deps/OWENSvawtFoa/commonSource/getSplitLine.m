function [data] = getSplitLine(fid,delim)

line = myfgetl(fid);
% Find where all of the delimiters are
delimiter_idx = find(line == delim);
delimiter_idx = [0.0,delimiter_idx,length(line)+1];
% Extract the data from the beginning to the last delimiter
data = zeros(length(delimiter_idx)-1,1);
for k = 2:length(delimiter_idx)
    data(k-1) = real(str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1)));
end

end
