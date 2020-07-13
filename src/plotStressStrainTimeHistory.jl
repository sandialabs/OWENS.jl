function plotStressStrainTimeHistory(time,strainHistory,elNum,gpNum,y,z,E,Gxy,Gxz)

    for i=2:length(time)
        strain = strainHistory(elNum,i-1);
        [eps_xx(i),gam_xy(i),gam_xz(i)] = calculateStrain(strain,gpNum,y,z);
    end
    
    figure(1);
    plot(time,eps_xx);
    grid on;
    xlabel('time (s)');
    ylabel('\epsilon_{xx}');
    
    figure(2);
    plot(time,gam_xy);
    grid on;
    xlabel('time (s)');
    ylabel('\gamma_{xy}');
    
    figure(3);
    plot(time,gam_xz);
    grid on;
    xlabel('time (s)');
    ylabel('\gamma_{xz}');
    
    figure(4);
    plot(time,E*eps_xx/1e6);
    grid on;
    xlabel('time (s)');
    ylabel('\sigma_{xx} (MPa)');
    
    figure(5);
    plot(time,Gxy*gam_xy/1e6);
    grid on;
    xlabel('time (s)');
    ylabel('\sigma_{xy} (MPa)');
    
    figure(6);
    plot(time,Gxz*gam_xz/1e6);
    grid on;
    xlabel('time (s)');
    ylabel('\sigma_{xz} (MPa)');
    
end