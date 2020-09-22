#######################################
### Verification/Regression Script ###
#######################################

# SETUP

tic()
verify_transient = true;
verify_modal = true;
plot_modal = true;
verify_flutter = true;

# RUN OWENS
test_owens(verify_transient,verify_modal,verify_flutter);

bmOwens = './input_files_test/_15mTower_transient_dvawt_c_2_lcdt';
[mesh,meshSeg] = readMeshVG([bmOwens '.mesh']); #read mesh file
tol1 = 1e-5;

#MODAL

if verify_modal
    OLD = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out');
    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out');

    FileInfo = dir(NEW);

    if (datetime - FileInfo.date) > duration(0,1,0)
        warning('Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.');
    end

    modalOutold = readResultsModalOut(OLD,mesh.numNodes);
    modalOutnew = readResultsModalOut(NEW,mesh.numNodes);



    fprintf('#s\n','MODAL TESTS ');
    total_mismatch = 0;
    for i = 1:2:20
        selectedMode = i;
        fprintf('\n#s','MODE: ');
        fprintf('#0.0f\n',i);

        # Get and unitize the selected modeshapes
        modeShapeold = modalOutold{selectedMode}.InPhaseShape;
        modalOutold{selectedMode}.InPhaseShape = table2struct(normalizemode(modeShapeold));
        modeShapeOutPhaseold = modalOutold{selectedMode}.OutPhaseShape;
        modalOutold{selectedMode}.OutPhaseShape = table2struct(normalizemode(modeShapeOutPhaseold));

        modeShapenew = modalOutnew{selectedMode}.InPhaseShape;
        modalOutnew{selectedMode}.InPhaseShape = table2struct(normalizemode(modeShapenew));
        modeShapeOutPhasenew = modalOutnew{selectedMode}.OutPhaseShape;
        modalOutnew{selectedMode}.OutPhaseShape = table2struct(normalizemode(modeShapeOutPhasenew));

        # Check for mismatch
        total_mismatch = total_mismatch + comparison(modalOutold{selectedMode},modalOutnew{selectedMode},tol1);

    end

    toc

    if total_mismatch == 0
        fprintf('#s\n\n','Modal Tests PASSED')
    else
        fprintf('#s\n\n','!!!Modal Tests FAILED!!!')
    end

end

if plot_modal
    disp('Plotting Modes')
    Ndof = 10;
    savePlot = true;

    n_comparisons = 4;
    outnameC = cell(1,n_comparisons);
    outnameC{1} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_CPP.out';
    outnameC{2} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION_EIGS.out';
    outnameC{3} = './input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out';
    outnameC{4} = './input_files_test/SORTED_EIG_VERIF_TEMP.out';

    for ii = 1:n_comparisons
        outname = outnameC{ii};
        disp(outname);
        if 1 # run to pause through plots
            for df = 1:2:Ndof
                viz([bmOwens '.mesh'],outname,df,10)
                set(gcf,'visible','off')
                if savePlot # save the plot
                    saveas(gcf,[outname(1:end-4) '_MODE' num2str(df) '.pdf'])
                    close gcf
                else # flip through the plots visually
                    pause
                end
            end
        else # generate all the plots and then view them
            for df = 1:Ndof
                df_act = Ndof-df+1;
                viz([bmOwens '.mesh'],outname,df_act,10)
                title(['DOF: ' num2str(df_act)])
            end
        end
    end
fprintf('#s\n\n','MODAL PLOTTING COMPLETE');

end


#FLUTTER
if verify_flutter

    fprintf('#s\n','FLUTTER TESTS ');

    freq_old = [0,8796.17977509167,0,8731.33585194320]';
    damp_old = [0,5.06872269198530e-09,0,6.95419907854427e-07]';
    old1.freq = freq_old;
    old1.damp = damp_old;

    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_FLUTTER.out');

    FileInfo = dir(NEW);

    if (datetime - FileInfo.date) > duration(0,1,0)
        warning('Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.');
    end

    new = load(NEW);
    freq_new = new(:,1);
    damp_new = new(:,2);
    new1.freq = freq_new;
    new1.damp = damp_new;

    total_mismatch = comparison(old1,new1,1e-3);

    if total_mismatch == 0
        fprintf('#s\n\n','Flutter Test PASSED')
    else
        fprintf('#s\n\n','!!!Flutter Test FAILED!!!')
    end

end

#TRANSIENT

if verify_transient

    n_t = 50;

    OLD = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff0.mat');
    NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.mat');

    FileInfo = dir(NEW);

    if (datetime - FileInfo.date) > duration(0,1,0)
        warning('Output was not generated, cannot compare stale output, a recent change must have prevented the output from being written or read in.');
    end


    # old = load(OLD);
    # new = load(NEW);

    old = loadOWENSmat(OLD,n_t);
    new = loadOWENSmat(NEW,n_t);
    fprintf('#s\n','TRANSIENT TESTS ');
    total_mismatch = comparison(old,new,tol1);

    toc

    if total_mismatch == 0
        fprintf('#s\n','Transient Tests PASSED')
    else
        fprintf('#s\n','!!!Transient Tests FAILED!!!')
    end

end

function [total_mismatch] = comparison(old,new,tol1)
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
                        msg = sprintf('#20s#1s#0.0i#2s#0.0i#6s#8e#6s#8e', subvarnames{j} ,' ', ii , ' : ' , jj ,' OLD: ' , sub_old_data(ii,jj) , ' NEW: ' , sub_new_data(ii,jj));
                        fprintf('#s\n',msg)
                        num_mismatch = num_mismatch + 1;
                    end
                end
            end
        end

    else#if ~(min(min(ismembertol(old_data,new_data,tol1))))
        for ii = 1:length(old_data(:,1))
            for jj = 1:length(old_data(1,:))
                if (abs(old_data(ii,jj)-new_data(ii,jj))>tol1)
                    msg = sprintf('#20s#1s#0.0f#2s#0.0f#6s#8e#6s#8e',varnames{i} ,' ', ii , ' : ' , jj ,' OLD: ' , old_data(ii,jj) , ' NEW: ' , new_data(ii,jj));
                    fprintf('#s\n',msg)
                    num_mismatch = num_mismatch + 1;
                end
            end
        end

    end

    if num_mismatch == 0
        fprintf('#s\n',['PASSED | ' , varnames{i}])
    else
        fprintf('#s\n',['FAILED | ', varnames{i}])
    end

    total_mismatch = total_mismatch + num_mismatch;

end

end

function [shape] = normalizemode(shape)
col_names = shape.Properties.VariableNames;
for i = 1:length(col_names)
    shape.(col_names{i}) = shape.(col_names{i})/max(abs(shape.(col_names{i})));
end

end
