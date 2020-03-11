function data = importCactusFile(filename,skiplines,row_len,col_len,delim)
fid = fopen(filename);
data = NaN(row_len,col_len); %TODO: dont make this hard coded
% skip header
for i = 1:skiplines
    myfgetl(fid);
end
j = 0;
while ~feof(fid)
    j = j+1;
    line = myfgetl(fid);
    
    % Find where all of the delimiters are
    delimiter_idx = [0.0,find(line == delim),length(line)+1];
    % Extract the data from the beginning to the last delimiter
    for k = 2:length(delimiter_idx)
        data(j,k-1) = str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1));
    end
end
fclose(fid);

data = data(1:j-1,:); %trim the excess off

end