function findStrainTimeExtrema(time,strainHistory,y,z)

    [numEl,~] = size(strainHistory);
    numGP = length(strainHistory(1).eps_xx_0);
    
    
    for i=2:length(time)
        index = 1;
        for j=1:numEl
            for k=1:numGP
            strain = strainHistory(j,i-1);
            [eps_xx(i,index),gam_xy(i,index),gam_xz(i,index)] = calculateStrain(strain,k,y,z);
            index = index + 1;
            end
        end
        
        [min_eps_xx(i),min_eps_xx_id(i)] = min(eps_xx(i,:));
        [max_eps_xx(i),max_eps_xx_id(i)] = max(eps_xx(i,:));
                
        [min_gam_xy(i),min_gam_xy_id(i)] = min(gam_xy(i,:));
        [max_gam_xy(i),max_gam_xy_id(i)] = max(gam_xy(i,:));
        
        [min_gam_xz(i),min_gam_xz_id(i)] = min(gam_xz(i,:));
        [max_gam_xz(i),max_gam_xz_id(i)] = max(gam_xz(i,:));
        
    end
    figure(1);
    subplot(2,1,1);
    plot(time,min_eps_xx,time,max_eps_xx);
    subplot(2,1,2);
    plot(time,min_eps_xx_id,time,max_eps_xx_id);
    
    figure(2);
    subplot(2,1,1);
    plot(time,min_gam_xy,time,max_gam_xy);
    subplot(2,1,2);
    plot(time,min_gam_xy_id,time,max_gam_xy_id);
    
    figure(3);
    subplot(2,1,1);
    plot(time,min_gam_xz,time,max_gam_xz);
    subplot(2,1,2);
    plot(time,min_gam_xz_id,time,max_gam_xz_id);
disp('test');
end