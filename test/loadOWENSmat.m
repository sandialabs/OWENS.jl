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
data1 = zeros(1,17);
for i = 0:n_t+1
    if i == 0
        line = myfgetl(fid);
    else
        data1 = getSplitLine(fid,',');

        out.t(i) = data1(1);
        out.aziHist(i) = data1(2);
        out.OmegaHist(i) = data1(3);
        out.OmegaDotHist(i) = data1(4);
        out.gbHist(i) = data1(5);
        out.gbDotHist(i) = data1(6);
        out.gbDotDotHist(i) = data1(7);
        out.FReactionHist(i,:) = data1(8:13);
        out.rigidDof(i) = data1(14);
        out.genTorque(i) = data1(15);
        out.genPower(i) = data1(16);
        out.torqueDriveShaft(i) = data1(17);

    end
end

fclose(fid);

% read uHist
fid = fopen([fname(1:end-4) '_uHist.txt'],'r');
data2 = zeros(1,492);
for i = 0:n_t+1
    if i == 0
        line = myfgetl(fid);
    else
        data2 = getSplitLine(fid,' ');
        out.uHist(:,i) = data2(2:end); % first is time
    end
end

fclose(fid);


% read strainHist

fid = fopen([fname(1:end-4) '_strainHist.txt'],'r');
data3 = zeros(1,4);
for i = 0:n_t
    if i == 0
        line = myfgetl(fid); % Skip Header
        line = myfgetl(fid);
    else

        for j = 1:length(out.strainHist(:,1))
            line = myfgetl(fid); % Skip Sub-Header
            line = myfgetl(fid);

            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_0(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_z(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).eps_xx_y(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xz_0(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xz_y(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xy_0(1,:) = data3(2:end);
            data3 = getSplitLine(fid, ' ');
            out.strainHist(j,i).gam_xy_z(1,:) = data3(2:end);
        end
    end
end

fclose(fid);
end
