% Based on Kattelans, Conservation of mass and momentum of the
% least-squares spectral collocation scheme for the Stokes problem

clear all
close all
clc


global N xi

N = 3;

xi = GLLnodes(N);

xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0.75 3 -1.5 -0.75 0 0.75 3 ...
            -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 ];
yglobal = [ -0.75 -0.75 -0.75 -0.75 -0.75 0 0 0 0 0.75 0.75 0.75 0.75 0.75 ...
            0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 ];
        
plot(xglobal,yglobal,'o','markerface','b')
hold on
axis equal
axis([-1.5 3 -0.75 0.75])


corners = [  1  2  7  6
             2 16 15  7
             2  3 17 16
             3  4 18 17
            18  4  8 19
             4  5  9  8
             6  7 11 10
             7 15 22 11
            22 21 12 11
            21 20 13 12
            19  8 13 20
             8  9 14 13 ];

for i=[1 6 7 12]
    
    X = transfinitemapping_v2(xglobal(corners(i,:)));
    Y = transfinitemapping_v2(yglobal(corners(i,:)));

    pcolor(reshape(X',N+1,N+1),reshape(Y',N+1,N+1),(i-1)/11*ones(N+1))
    
end



tblr = zeros(8,4);
tblr(1,4) = 1; tblr(2,1) = 1; tblr(3,1)  = 1; tblr(4,3)  = 1;
tblr(5,4) = 1; tblr(6,2) = 1; tblr(7,2) = 1; tblr(8,3) = 1;

theta = [ -3/4*pi     -pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
               pi  3/4*pi
           3/4*pi    pi/2
             pi/2    pi/4
                0    pi/4 ];

R = 0.5;

k=0;
for i=[2:5 8:11]

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];
    [X,Y] = transfinitemapping_cilinder_v2(cor,R,theta(k,:),tblr(k,:));
    pcolor(X,Y,(i-1)/11*ones(N+1))

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cilinder

Rc = R;
xc = [linspace(-1/2,1/2,100) linspace(1/2,-1/2,100)]; % coordinates of cilinder
yc = sqrt(Rc^2-xc.^2); yc(101:200) = -yc(101:200);

plot(xc,yc,'k','linewidth',4)







