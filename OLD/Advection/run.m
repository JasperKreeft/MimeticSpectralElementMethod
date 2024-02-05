% Curves problem
clear all
close all
clc

N = 3;

x0 = linspace(0,1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

dT = 1;

for i=1:N
    
    [t,y] = ode45(@f1,[0 dT],x0(i));

    plot(y,t,'-x'); hold on
end

grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

dT = 2*pi;

y0 = 0.0;

for i=1:N
   
[t,xy] = ODE45(@f2,[0 dT],[x0(i) y0]);

x = xy(:,1);
y = xy(:,2);

plot(x,y,'.-'); hold on

end

grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

ds = 10*pi;

t0 = 0.0;
x0 = 0.5;
y0 = 0.0;

[s,txy] = ODE45(@f3,[0 ds],[t0 x0 y0]);

t = txy(:,1);
x = txy(:,2);
y = txy(:,3);

plot3(x,y,t,'-')
grid on