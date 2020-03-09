function result = matc(sf,tf,tol1,showdiff)
    % Compare two mat files
    % matc matfile1.mat matfile2.mat
    
    % Copyright 2004, Dr. Xianyao Chen, xianyaochen@yahoo.com.cn
    % Last Modified by Xianyao Chen, 08-June-2004
    % From https://www.mathworks.com/matlabcentral/fileexchange/5160-compare-two-mat-files
    
    N_not_equal = 0;
    
%     if (nargin ~= 2)
%         return
%     end
%     if exist(sf) ~=2
%         msg = sprintf('%s%s','Not valid mat file: ',sf);
%         %     fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
%         fprintf('%s',msg)
%         return
%     end
%     if exist(tf) ~=2
%         msg = sprintf('%s%s','Not valid mat file: ',tf);
%         %     fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
%         fprintf('%s',msg)
%         return
%     end
    save #~tmp1.mat -mat
    clear
    load #~tmp1 sf
    
    load(sf)
    clear sf
    
    ws = whos;
    save #~tmp2.mat -mat
    clear
    
    load #~tmp1
    load #~tmp2
    delete('#~tmp2.mat');
    
    m1 = size(ws,1);
    
    % warning off
    % nPos = 0;
    for k = 1:m1
        var = ws(k).name;
        as = load(sf,var);
        at = load(tf,var);
        atn = fieldnames(at);
        if k == 1
            msg = sprintf('%s%s','Source File: ',sf,' --- Target File: ',tf');
            fprintf('%s\n',msg)
            msg = sprintf('%50s','--------------------------------------------------');
            fprintf('%s\n',msg)
            %         msg = sprintf('%s','----------------------------------------------------------------------------------------------------------------------------');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
        end
        
        if isempty(atn)
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' not exist');
            fprintf('%s\n',msg)
            %         msg = sprintf('%12s%12s%20s%20s%20s%20s%s','Variable :: ',var,'    |  Source File: ',sf,' ---  Target File: ',tf,' :: not exist');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            continue
        end
        aas = struct2cell(as);
        aat = struct2cell(at);
        
        % Determine if equal
        if isnumeric(aas{1}) && isnumeric (aat{1})
            equal_logic = (min(min(ismembertol(aas{1},aat{1},tol1))));
        else
            equal_logic = isequal(aas{1},aat{1});
        end
        
        if ~equal_logic
            N_not_equal = N_not_equal+1;
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' not equal');
            fprintf('%s\n',msg)
            %         msg = sprintf('%12s%12s%20s%20s%20s%20s%s','Variable :: ',var,'    |  Source File: ',sf,' ---  Target File: ',tf,' :: not equal');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            
            if showdiff
                
                content1 = aas{1};
                content2 = aat{1};
                
                for ii = 1:length(content1(:,1))
                    for jj = 1:length(content1(1,:))
                        if ~isequal(content1,content2) && (isnumeric(content1) && isnumeric (content2))
                            msg = sprintf('%12s%i%2s%i%6s%8e%6s%8e','     Element ' , ii , ':' , jj ,' SRC: ' , content1(ii,jj) , ' TRG: ' , content2(ii,jj));
                            fprintf('%s\n',msg)
                        end
                    end
                end
            end
            
            continue
        else
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' equal');
            fprintf('%s\n',msg)
            %         msg = sprintf('%12s%12s%20s%20s%20s%20s%s','Variable :: ',var,'    |  Source File: ',sf,' ---  Target File: ',tf,' :: Exactly Equal');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            continue
        end
    end
    msg = sprintf('%50s','--------------------------------------------------');
    fprintf('%s\n\n',msg)
    
    % msg = sprintf('%s','----------------------------------------------------------------------------------------------------------------------------');
    % fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
    % fprintf('\n')
    
    
    clear
    load #~tmp1 tf
    
    load(tf)
    clear tf
    
    ws = whos;
    save #~tmp2.mat -mat
    clear
    
    load #~tmp1
    load #~tmp2
    delete('#~tmp2.mat');
    
    m1 = size(ws,1);
    % warning off
    
    for k = 1:m1
        var = ws(k).name;
        as = load(tf,var);
        at = load(sf,var);
        atn = fieldnames(at);
        if k == 1
            msg = sprintf('%s%s','Source File: ',tf,' --- Target File: ',sf');
            fprintf('%s\n',msg)
            msg = sprintf('%50s','--------------------------------------------------');
            fprintf('%s\n',msg)
            %fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %fprintf('\n')
        end
        
        if isempty(atn)
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' not exist');
            fprintf('%s\n',msg)
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            continue
        end
        aas = struct2cell(as);
        aat = struct2cell(at);
        
        % Determine if equal
        if isnumeric(aas{1}) && isnumeric (aat{1})
            equal_logic = (min(min(ismembertol(aas{1},aat{1},tol1))));
        else
            equal_logic = isequal(aas{1},aat{1});
        end
        
        if ~equal_logic
            N_not_equal = N_not_equal+1;
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' not equal');
            fprintf('%s\n',msg)
            %         msg = sprintf('%12s%12s%20s%20s%20s%20s%s','Variable :: ',var,'    |  Source File: ',tf,' ---  Target File: ',sf,' :: not equal');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            
            if showdiff
                content1 = aat{1};
                content2 = aas{1};
                
                for ii = 1:length(content1(:,1))
                    for jj = 1:length(content1(1,:))
                        if ~isequal(content1,content2) && (isnumeric(content1) && isnumeric (content2))
                            msg = sprintf('%12s%i%2s%i%6s%8e%6s%8e','     Element ' , ii , ':' , jj ,' SRC: ' , content1(ii,jj) , ' TRG: ' , content2(ii,jj));
                            fprintf('%s\n',msg)
                        end
                    end
                end
            end
            
            continue
        else
            msg = sprintf('%12s%12s%6s%s','Variable :: ',var,'  ||  ',' equal');
            fprintf('%s\n',msg)
            %         msg = sprintf('%12s%12s%20s%20s%20s%20s%s','Variable :: ',var,'    |  Source File: ',tf,' ---  Target File: ',sf,' :: Exactly Equal');
            %         fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
            %         fprintf('\n')
            continue
        end
    end
    % warning on
    
    msg = sprintf('%50s','--------------------------------------------------');
    fprintf('%s\n',msg)
    % fprintf( [ repmat( '\b', 1, nPos ), '%s' ], msg),
    % fprintf('\n')
    
    if N_not_equal==0
        msg = sprintf('%12s','Tests Passed');
        fprintf('%s\n',msg)
    else
        msg = sprintf('%12s','!!!TESTS FAILED!!!');
        fprintf('%s\n',msg)
        error('!!!TESTS FAILED Something Changed!!!');
    end
    
    
    delete('#~tmp1.mat')