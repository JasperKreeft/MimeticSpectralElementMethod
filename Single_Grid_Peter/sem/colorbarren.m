clear all
close all
clc

x = 1:10;
y = 1:2;

for i=1:length(x)
    for j=1:length(y)
        z(j,i)=i;
    end
end

subplot(2,1,1)
imagesc(x,y,z)


subplot(2,1,2)
pcolor(x,y,z)
shading interp