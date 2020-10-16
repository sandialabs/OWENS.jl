% C++ compatable workaround for fget1 from https://www.mathworks.com/matlabcentral/answers/461159-read-text-file-line-by-line-in-deployed-application
function line  = myfgetl(fid)
line = '';
c = fread(fid,1,'char=>char');
while ~feof(fid) && ~strcmp(c, sprintf('\n'))
    if ~strcmp(c,sprintf('\r'))
        line = [line,c];
    end
    c = fread(fid,1,'char=>char');
end
end
