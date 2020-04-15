tic()
verify_transient = false;
verify_modal = false;
plot_modal = true;
verify_flutter = false;

test_owens(verify_transient,verify_modal,verify_flutter);

if verify_modal
    OLD = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out');
    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out');
    
    old = fread(fopen(OLD));
    new = fread(fopen(NEW));
    
    
    
    if ~isequal(old,new)
        fprintf('%s\n','!!!MODAL TESTS FAILED!!!')
    else
        fprintf('%s\n','MODAL TESTS PASSED')
    end
end

if plot_modal
    disp('Plotting Modes')
    Ndof = 10;
    savePlot = true;
    
    bmOwens = './input_files_test/_15mTower_transient_dvawt_c_2_lcdt';
    n_comparisons = 4;
    outnameC = cell(1,n_comparisons);
    outnameC{1} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_CPP.out';
    outnameC{2} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION_EIGS.out';
    outnameC{3} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out';
    outnameC{4} = './input_files_test/SORTED_EIG_VERIF_TEMP.out';
    
    for ii = 1:n_comparisons
        outname = outnameC{ii};
        disp(outname);
        if 1 % run to pause through plots
            for df = 1:2:Ndof
                viz([bmOwens '.mesh'],outname,df,10)
                set(gcf,'visible','off')
                if savePlot % save the plot
                    saveas(gcf,[outname(1:end-4) '_MODE' num2str(df) '.pdf'])
                    close gcf
                else % flip through the plots visually
                    pause
                end
            end
        else % generate all the plots and then view them
            for df = 1:Ndof
                df_act = Ndof-df+1;
                viz([bmOwens '.mesh'],outname,df_act,10)
                title(['DOF: ' num2str(df_act)])
            end
        end
    end
    
    
end

if verify_transient
    
    n_t = 50;
    
    tol1 = 1e-5;
    OLD = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff0.mat');
    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.mat');
    
    % FileInfo = dir(NEW);
    
    % if (datetime - FileInfo.date) > duration(0,1,0)
    %     error('Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.');
    % end
    
    
    % old = load(OLD);
    % new = load(NEW);
    
    old = loadOWENSmat(OLD,n_t);
    new = loadOWENSmat(NEW,n_t);
    
    total_mismatch = 0;
    
    varnames = fieldnames(old);
    for i = 1:length(varnames)
        
        num_mismatch = 0;
        old_data = old.(varnames{i});
        new_data = new.(varnames{i});
        
        if ~isnumeric(old_data)
            subvarnames = fieldnames(old_data);
            for j = 1:length(subvarnames)
                
                sub_old_data = old_data.(subvarnames{j});
                sub_new_data = new_data.(subvarnames{j});
                for ii = 1:length(sub_old_data(:,1))
                    for jj = 1:length(sub_old_data(1,:))
                        if (abs(sub_old_data(ii,jj)-sub_new_data(ii,jj))>tol1)
                            msg = sprintf('%20s%i%2s%i%6s%8e%6s%8e', subvarnames{j}  , ii , ':' , jj ,' OLD: ' , sub_old_data(ii,jj) , ' NEW: ' , sub_new_data(ii,jj));
                            fprintf('%s\n',msg)
                            num_mismatch = num_mismatch + 1;
                        end
                    end
                end
            end
            
        else%if ~(min(min(ismembertol(old_data,new_data,tol1))))
            for ii = 1:length(old_data(:,1))
                for jj = 1:length(old_data(1,:))
                    if (abs(old_data(ii,jj)-new_data(ii,jj))>tol1)
                        msg = sprintf('%20s%f%2s%f%6s%8e%6s%8e',varnames{i} , ii , ':' , jj ,' OLD: ' , old_data(ii,jj) , ' NEW: ' , new_data(ii,jj));
                        fprintf('%s\n',msg)
                        num_mismatch = num_mismatch + 1;
                    end
                end
            end
            
        end
        
        if num_mismatch == 0
            fprintf('%s\n',['PASSED | ' , varnames{i}])
        else
            fprintf('%s\n',['FAILED | ', varnames{i}])
        end
        
        total_mismatch = total_mismatch + num_mismatch;
        
    end
    
    toc
    
    if total_mismatch == 0
        fprintf('%s\n','TRANSIENT TESTS PASSED')
    else
        fprintf('%s\n','!!!TRANSIENT TESTS FAILED!!!')
    end
end