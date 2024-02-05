% clear all
% close all
% clc

phi = reshape(phi,1,[]);
x = reshape(XEGEG,1,[]);
y = reshape(YEGEG,1,[]);

M = length(x);

% x = linspace(0,1,M);
% y = x;
% 
% phi = linspace(0,1,M);

phi_min = -1;%min(phi);
phi_max = 1;%max(phi);


clim = [phi_min phi_max];

colormap;

colormap_phi = zeros(M,3);
for i=1:M
    colorindicator = (phi(i)-phi_min)/(phi_max-phi_min)*64;
    
    if colorindicator<=8
        colormap_phi(i,:) = [ 0 0 1/2+colorindicator/16];
    elseif colorindicator<=24
        colormap_phi(i,:) = [0 (colorindicator-8)/16 1];
    elseif colorindicator<=40
        colormap_phi(i,:) = [(colorindicator-24)/16 1 (40-colorindicator)/16];
    elseif colorindicator<=56
        colormap_phi(i,:) = [1 (56-colorindicator)/16 0];
    elseif colorindicator<=64
        colormap_phi(i,:) = [1/2+(64-colorindicator)/16 0 0];
    end
    
end

scatter3(x,y,phi,100*ones(M,1),colormap_phi,'filled')
colorbar('yticklabel',roundn(linspace(phi_min,phi_max,6),-2))