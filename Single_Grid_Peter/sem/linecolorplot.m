function []=linecolorplot(q,x,y,qmin,qmax)

% q = reshape(q,1,[]);
% x = reshape(x,1,[]);
% y = reshape(y,1,[]);

if nargin<=3
    qmin = min(q);
    qmax = max(q);
end

M = length(x);

clim = [qmin qmax];

colormap;

colormap_q = zeros(M,3);
for i=1:M
    colorindicator = (q(i)-qmin)/(qmax-qmin)*64;
    
    if colorindicator<=8
        colormap_q(i,:) = [ 0 0 1/2+colorindicator/16];
    elseif colorindicator<=24
        colormap_q(i,:) = [0 (colorindicator-8)/16 1];
    elseif colorindicator<=40
        colormap_q(i,:) = [(colorindicator-24)/16 1 (40-colorindicator)/16];
    elseif colorindicator<=56
        colormap_q(i,:) = [1 (56-colorindicator)/16 0];
    elseif colorindicator<=64
        colormap_q(i,:) = [1/2+(64-colorindicator)/16 0 0];
    end

end

for i=1:M
plot3(x(i,:),y(i,:),q(i,:),'color',colormap_q(i,:),'linewidth',4)
hold on
end
colorbar('yticklabel',roundn(linspace(qmin,qmax,6),-2))