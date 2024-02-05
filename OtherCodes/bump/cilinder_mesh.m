cilinder;

N = 8;

% xi = linspace(-1,1,N+1)';
xi = GLLnodes(N)';
eta = xi';

Xi  = xi*ones(1,N+1);
Eta = Xi';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% block links

x = [Xol xv(3,1:3)']';
y = [Yol yv(3,1:3)']';

siz = size(x)-[1 1];

for i=1:siz(1)
    for j=1:siz(2)
        kl=i+(j-1)*siz(1);

        X = transfinitemapping(x,xi,eta,i,j);

        Y = transfinitemapping(y,xi,eta,i,j);
        
        figure(3); subplot(2,1,2)
        surf(X,Y,kl*ones(N+1))
        view([0 0 1])
        hold on
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% block rechts

x = [Xor xv(3,7:-1:5)']';
y = [Yor yv(3,7:-1:5)']';

siz = size(x)-[1 1];

for i=1:siz(1)
    for j=1:siz(2)
        kr=i+(j-1)*siz(1)+kl;

        X = transfinitemapping(x,xi,eta,i,j);

        Y = transfinitemapping(y,xi,eta,i,j);
        
        figure(3); subplot(2,1,2)
        surf(X,Y,kr*ones(N+1))
        view([0 0 1])
        hold on
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm     = sqrt(xv.^2+yv.^2)';
thetam = fliplr(atand(yv./xv)+(xv<=0)*180)';



siz = size(rm)-[1 1];

for i=1:siz(1)
    for j=1:siz(2)
        k=i+(j-1)*siz(1)+kr;

        Tm = transfinitemapping(thetam,xi,eta,i,j);

%         Tm = transfinitemapping_cilinder(rm(i:i+1,j:j+1),thetam(i:i+1,j:j+1),xi,eta,[0 0 0 0]);
        
        if j==1
            tblr = [0 0 0 0];
        elseif j==2
            tblr = [1 0 0 0];%[1 1 0 0];
        end
%         Rm = coordinate(rm,xi,eta,i,j);
        Rm = transfinitemapping_cilinder(rm(i:i+1,j:j+1),thetam(i:i+1,j:j+1),xi,eta,tblr);
     
        Xm = Rm.*cosd(Tm);
        Ym = Rm.*sind(Tm);

        figure(3)
        subplot(2,1,1)
        surf(Tm,Rm,k*ones(N+1))
        axis([0 180 0 2.5])
        view([0 0 1])
        hold on

        subplot(2,1,2)
        surf(Xm,Ym,k*ones(N+1))
        axis equal
        axis([-4 4 0 2.5])
        view([0 0 1])
        hold on
    end
end


