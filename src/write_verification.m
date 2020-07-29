function write_verification(model,t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,rigidDof,genTorque,genPower,torqueDriveShaft,uHist,strainHist)


    fprintf('#s\n','Writing Verification File')


    fid = fopen([model.outFilename(1:end-3) 'txt'],'w+'); %open/create new for writing and discard existing data

    for i = 0:length(t)
        if i == 0
            line = ['t,aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist1,FReactionHist2,FReactionHist3,FReactionHist4,FReactionHist5,FReactionHist6,rigidDof,genTorque,genPower,torqueDriveShaft' sprintf('\n')'];
        else
            line = [sprintf('%1.12e', t(i)), ',',...
            sprintf('%1.12e', aziHist(i)), ',',...
            sprintf('%1.12e', OmegaHist(i)), ',',...
            sprintf('%1.12e', OmegaDotHist(i)), ',',...
            sprintf('%1.12e', gbHist(i)), ',',...
            sprintf('%1.12e', gbDotHist(i)), ',',...
            sprintf('%1.12e', gbDotDotHist(i)), ',',...
            sprintf('%1.12e', FReactionHist(i,1)), ',',...
            sprintf('%1.12e', FReactionHist(i,2)), ',',...
            sprintf('%1.12e', FReactionHist(i,3)), ',',...
            sprintf('%1.12e', FReactionHist(i,4)), ',',...
            sprintf('%1.12e', FReactionHist(i,5)), ',',...
            sprintf('%1.12e', FReactionHist(i,6)), ',',...
            sprintf('%1.12e', rigidDof(i)), ',',...
            sprintf('%1.12e', genTorque(i)), ',',...
            sprintf('%1.12e', genPower(i)), ',',...
            sprintf('%1.12e', torqueDriveShaft(i)), sprintf('\n') ];
        end
        fwrite(fid,line,'char');
    end

    fclose(fid);

    % Save uHist
    fid = fopen([model.outFilename(1:end-4) '_uHist.txt'],'w+'); %open/create new for writing and discard existing data

    for i = 0:length(t)
        if i == 0
            line = ['t uHist' sprintf('\n')'];
            fwrite(fid,line,'char');
        else
            for ii = 1:length(uHist(:,1))-1
                line = [sprintf('%1.12e',t(i)), ' ', sprintf('%1.12e',uHist(ii,i)') ];
                fwrite(fid,line,'char');
            end
            line = [sprintf('%1.12e',t(i)), ' ', sprintf('%1.12e',uHist(end,i)'), sprintf('\n') ];
            fwrite(fid,line,'char');
        end

    end

    fclose(fid);


    % Save strainHist
    fid = fopen([model.outFilename(1:end-4) '_strainHist.txt'],'w+'); %open/create new for writing and discard existing data

    for i = 0:length(t)-1
        if i == 0
            line = ['n_t n_elem' sprintf('\n')];
            fwrite(fid,line,'char');

            line = [sprintf('%1.12e',length(t)) ' ' sprintf('%1.12e',length(strainHist(:,1))) sprintf('\n')];
            fwrite(fid,line,'char');
        else

            for j = 1:length(strainHist(:,1))
                line = ['t elem', sprintf('\n')];
                fwrite(fid,line,'char');
                line = [sprintf('%1.12e',t(i)), ' ', sprintf('%1.12e',j), sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['eps_xx_0 ' sprintf('%1.12e',strainHist{j,i}.eps_xx_0(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['eps_xx_0 ' sprintf('%1.12e',strainHist{j,i}.eps_xx_0(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['eps_xx_z ' sprintf('%1.12e',strainHist{j,i}.eps_xx_z(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['eps_xx_z ' sprintf('%1.12e',strainHist{j,i}.eps_xx_z(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['eps_xx_y ' sprintf('%1.12e',strainHist{j,i}.eps_xx_y(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['eps_xx_y ' sprintf('%1.12e',strainHist{j,i}.eps_xx_y(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['gam_xz_0 ' sprintf('%1.12e',strainHist{j,i}.gam_xz_0(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['gam_xz_0 ' sprintf('%1.12e',strainHist{j,i}.gam_xz_0(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['gam_xz_y ' sprintf('%1.12e',strainHist{j,i}.gam_xz_y(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['gam_xz_y ' sprintf('%1.12e',strainHist{j,i}.gam_xz_y(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['gam_xy_0 ' sprintf('%1.12e',strainHist{j,i}.gam_xy_0(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['gam_xy_0 ' sprintf('%1.12e',strainHist{j,i}.gam_xy_0(4)) sprintf('\n')];
                fwrite(fid,line,'char');

                for jj = 1:3
                    line = ['gam_xy_z ' sprintf('%1.12e',strainHist{j,i}.gam_xy_z(jj))];
                    fwrite(fid,line,'char');
                end
                line = ['gam_xy_z ' sprintf('%1.12e',strainHist{j,i}.gam_xy_z(4)) sprintf('\n')];
                fwrite(fid,line,'char');

            end
        end
    end

    fclose(fid);
end
