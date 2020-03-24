function [data] = getSplitLine(fid,delim)

line = myfgetl(fid);

% Find where all of the delimiters are
delimiter_idx = [0.0,find(line == delim),length(line)+1];

% % Reduce index, getting rid of duplicate delimiters if spaces
% if delim == ' '
%     use_idx = true(1,length(delimiter_idx));
%     for i = 1:length(delimiter_idx)
%         for j = 1:length(delimiter_idx)-i
%             if delimiter_idx(i)+j == delimiter_idx(i+j)
%                 use_idx(i) = false;
%             else
%                 break
%             end
%         end
%     end
%     delimiter_idx = delimiter_idx(use_idx);
% end


% Extract the data from the beginning to the last delimiter
data = zeros(length(delimiter_idx)-1,1);
for k = 2:length(delimiter_idx)
    data(k-1) = real(str2double(line(delimiter_idx(k-1)+1:delimiter_idx(k)-1)));
end

end
