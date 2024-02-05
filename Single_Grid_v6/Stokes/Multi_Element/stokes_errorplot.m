close all

axisXY = [ 1e-2 4e-1 1e-10 1];

kleur = 'brgc';
symb = 'o^sv';

for i=1:4

str1 = ['-' symb(i) kleur(i)];
str2 = kleur(i);
str3 = ['--' kleur(i)];

filename = ['SIAM_N' num2str(2*i) '.mat'];

load(filename)

ax = 1./(HconvRange);

subplot(2,2,2)
loglog(ax,Hd_w,str1,'markerface',str2)
hold on
% loglog(ax,L2_w,str3,'markerface',str2)
axis(axisXY)
title('w error Stokes')
set(gca,'ytick',10.^(-10:2:0))
legend('N=2','N=4','N=6','N=8',4)
set(legend,'box','off')
subplot(2,2,3)
loglog(ax,Hd_u,str1,'markerface',str2)
hold on
% loglog(ax,L2_U,str3,'markerface',str2)
axis(axisXY)
title('u error Stokes')
set(gca,'ytick',10.^(-10:2:0))
legend('N=2','N=4','N=6','N=8',4)
set(legend,'box','off')
subplot(2,2,4)
loglog(ax,L2_p,str1,'markerface',str2)
hold on
axis(axisXY)
title('p error Stokes')
set(gca,'ytick',10.^(-10:2:0))
legend('N=2','N=4','N=6','N=8',4)
set(legend,'box','off')
end