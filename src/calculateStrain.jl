function [eps_xx,gam_xy,gam_xz] = calculateStrain(strain,gpNum,y,z)
	
	eps_xx = strain.eps_xx_0(gpNum) + strain.eps_xx_y(gpNum)*y + strain.eps_xx_z(gpNum)*z;
    gam_xy = strain.gam_xy_0(gpNum) + strain.gam_xy_z(gpNum)*z;
    gam_xz = strain.gam_xz_0(gpNum) + strain.gam_xz_y(gpNum)*y;
		
end