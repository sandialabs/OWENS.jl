function out = loadOWENSmat(fname,n_t)

%readfile
out = struct('uHist',zeros(492,n_t+1),...
't',zeros(1,n_t+1),...
'aziHist',zeros(1,n_t+1),...
'OmegaHist',zeros(1,n_t+1),...
'OmegaDotHist',zeros(1,n_t+1),...
'gbHist',zeros(1,n_t+1),...
'gbDotHist',zeros(1,n_t+1),...
'gbDotDotHist',zeros(1,n_t+1),...
'FReactionHist',zeros(n_t+1,6),...
'rigidDof',zeros(1,n_t+1),...
'genTorque',zeros(1,n_t+1),...
'genPower',zeros(1,n_t+1),...
'torqueDriveShaft',zeros(1,n_t+1));

single_strainHist = struct('eps_xx_0',zeros(1,4),...
'eps_xx_z',zeros(1,4),...
'eps_xx_y',zeros(1,4),...
'gam_xz_0',zeros(1,4),...
'gam_xz_y',zeros(1,4),...
'gam_xy_0',zeros(1,4),...
'gam_xy_z',zeros(1,4));

out.strainHist = repmat(single_strainHist,75,n_t);

fid = fopen([fname(1:end-3) 'txt'],'r');

for i = 0:n_t+1
    if i == 0
        line = myfgetl(fid);
    else
        data = getSplitLine(fid,',');

        out.t(i) = data(1);
        out.aziHist(i) = data(2);
        out.OmegaHist(i) = data(3);
        out.OmegaDotHist(i) = data(4);
        out.gbHist(i) = data(5);
        out.gbDotHist(i) = data(6);
        out.gbDotDotHist(i) = data(7);
        out.FReactionHist(i,:) = data(8:13);
        out.rigidDof(i) = data(14);
        out.genTorque(i) = data(15);
        out.genPower(i) = data(16);
        out.torqueDriveShaft(i) = data(17);

    end
end

fclose(fid);

% read uHist
fid = fopen([fname(1:end-4) '_uHist.txt'],'r');

for i = 0:n_t+1
    if i == 0
        line = myfgetl(fid);
    else
        data = getSplitLine(fid,' ');
        out.uHist(:,i) = data(2:end); % first is time
    end
end

fclose(fid);


% read strainHist

fid = fopen([fname(1:end-4) '_strainHist.txt'],'r');

for i = 0:n_t
    if i == 0
        line = myfgetl(fid); % Skip Header
        line = myfgetl(fid);
    else

        for j = 1:length(out.strainHist(:,1))
            line = myfgetl(fid); % Skip Sub-Header
            line = myfgetl(fid);

            eps_xx_0 = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_0(1,:) = eps_xx_0(2:end);
            eps_xx_z = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_z(1,:) = eps_xx_z(2:end);
            eps_xx_y = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_y(1,:) = eps_xx_y(2:end);
            gam_xz_0 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xz_0(1,:) = gam_xz_0(2:end);
            gam_xz_y = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xz_y(1,:) = gam_xz_y(2:end);
            gam_xy_0 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xy_0(1,:) = gam_xy_0(2:end);
            gam_xy_z = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xy_z(1,:) = gam_xy_z(2:end);
        end
    end
end

fclose(fid);
end
