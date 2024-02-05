
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load centerlines_velocity_Reference.mat

figure
set(gcf,'Position',[1 46 280 925])
subplot(4,1,1)
plot(xH,uH,'r')
hold on
axis([-1 1 -.3 0])
subplot(4,1,2)
plot(xH,vH,'r')
hold on
axis([-1 1 -.2 .2])
subplot(4,1,3)
plot(uV,yV,'r')
hold on
axis([-.5 1 -1 1])
subplot(4,1,4)
plot(vV,yV,'r')
hold on
axis([-.1 .1 -1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('numRows','var')
    numRows    = 1;
    numColumns = 1;
end

s = size(Xp,1);

xH = zeros(1,numColumns*s);
uH = zeros(1,numColumns*s);
vH = zeros(1,numColumns*s);
ind1 = numColumns*floor(numRows/2);
ind2 = floor((s+1)/2);
ind3 = ceil((s+1)/2);
for i=1:numColumns
    ind = (i-1)*s+(1:s);
    xH(1,ind) = (Xp(:,ind2,i+ind1)+Xp(:,ind3,i+ind1))/2;
    uH(1,ind) = (uu(:,ind2,i+ind1)+uu(:,ind3,i+ind1))/2;
    vH(1,ind) = (vv(:,ind2,i+ind1)+vv(:,ind3,i+ind1))/2;
end

yV = zeros(numRows*s,1);
uV = zeros(numRows*s,1);
vV = zeros(numRows*s,1);
ind1 = ceil(numColumns/2);
for i=1:numRows
    ind = (i-1)*s+(1:s);
    yV(ind,1) = (Yp(ind2,:,(i-1)*numColumns+ind1)+Yp(ind3,:,(i-1)*numColumns+ind1))/2;
    uV(ind,1) = (uu(ind2,:,(i-1)*numColumns+ind1)+uu(ind3,:,(i-1)*numColumns+ind1))/2;
    vV(ind,1) = (vv(ind2,:,(i-1)*numColumns+ind1)+vv(ind3,:,(i-1)*numColumns+ind1))/2;
end


subplot(4,1,1)
plot(xH,uH)
subplot(4,1,2)
plot(xH,vH)
subplot(4,1,3)
plot(uV,yV)
subplot(4,1,4)
plot(vV,yV)
