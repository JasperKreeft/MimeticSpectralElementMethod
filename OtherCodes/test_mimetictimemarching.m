clear all
clf%ose all
clc

T = 10;
t = 0;
j=1;
dt = 0.05;

q(1) = 10;
while t<T
    t = t+dt;
    j = j+1;
    q(j) = mimetictimemarching(1,1,dt,q(j-1),3);
end

plot(linspace(0,t,j),q,'-o')
hold on
t=linspace(0,10,100);
plot(t,q(1)*exp(-t),'r')