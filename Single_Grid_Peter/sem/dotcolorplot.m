function []=dotcolorplot(phi,x,y,phimin,phimax)

phi = reshape(phi,1,[]);
x = reshape(x,1,[]);
y = reshape(y,1,[]);

if nargin<=3
    phimin = min(phi);
    phimax = max(phi);
end

M = length(x);

clim = [phimin phimax];

colormap;

colormap_phi = zeros(M,3);
for i=1:M
    colorindicator = (phi(i)-phimin)/(phimax-phimin)*64;
    
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
colorbar('yticklabel',roundn(linspace(phimin,phimax,6),-2))