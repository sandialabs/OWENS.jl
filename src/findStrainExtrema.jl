function [min_eps_xx_0,min_gam_xy_0,min_gam_xz_0,max_eps_xx_0,max_gam_xy_0,max_gam_xz_0] = findStrainExtrema(strain)

#unpack strain components
numEl = length(strain);
[r,numGPPerEl] = size(strain(1).eps_xx_0);
eps_xx_0 = zeros(numEl*numGPPerEl,1);
gam_xz_0 = zeros(numEl*numGPPerEl,1);
gam_xy_0 = zeros(numEl*numGPPerEl,1);

index = 1;
for i=1:numEl
    for j=1:numGPPerEl
        eps_xx_0(index,1) = strain(i).eps_xx_0(j);
        eps_xx_0(index,2) = i;
        eps_xx_0(index,3) = j;
        gam_xz_0(index,1) = strain(i).gam_xz_0(j);
        gam_xz_0(index,2) = i;
        gam_xz_0(index,3) = j;
        gam_xy_0(index,1) = strain(i).gam_xz_0(j);
        gam_xy_0(index,2) = i;
        gam_xy_0(index,3) = j;
        index = index + 1;
    end
end

[min_eps_xx_0,map1] = sort(eps_xx_0(:,1),'ascend');
min_eps_xx_0 = [min_eps_xx_0, eps_xx_0(map1,2),eps_xx_0(map1,3)];

[min_gam_xy_0,map2] = sort(gam_xy_0(:,1),'ascend');
min_gam_xy_0 = [min_gam_xy_0, gam_xy_0(map2,2),gam_xy_0(map2,3)];

[min_gam_xz_0,map3] = sort(gam_xz_0(:,1),'ascend');
min_gam_xz_0 = [min_gam_xz_0, gam_xz_0(map3,2),gam_xz_0(map3,3)];


[max_eps_xx_0,map4] = sort(eps_xx_0(:,1),'descend');
max_eps_xx_0 = [max_eps_xx_0, eps_xx_0(map4,2),eps_xx_0(map4,3)];

[max_gam_xy_0,map5] = sort(gam_xy_0(:,1),'descend');
max_gam_xy_0 = [max_gam_xy_0, gam_xy_0(map5,2),gam_xy_0(map5,3)];

[max_gam_xz_0,map6] = sort(gam_xz_0(:,1),'descend');
max_gam_xz_0 = [max_gam_xz_0, gam_xz_0(map6,2),gam_xz_0(map6,3)];


end
