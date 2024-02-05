% Triangle element build up from 3 quadrilateral elements

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N xi

N = 14;

xi = GLLnodes(N);

a = sqrt(3)/3;

xglobal = [ -1  0  1 1/2 0  -1/2 0 ];
yglobal = [ -a -a -a a/2 2*a a/2 0 ];


corners = [  1 2 7 6
             2 3 4 7
             7 4 5 6 ];

element1 = 1:(N+1)^2;
el2 = [ (N+1):(N+1):(N+1)^2
  reshape((N+1)^2+(1:N*(N+1)),N,N+1) ];
element2 = reshape(el2,1,(N+1)^2);

el3 = [ (N+1)^2-(0:N)
        [ (N+1)^2+N^2+(1:N)' reshape((N+1)^2+N*(N+1)+(1:N*N)',N,N)] ];
element3  = reshape(el3,1,(N+1)^2);
         
connectivity = [ element1
                 element2
                 element3 ];


% connectivity = [  1  2  3  4  5  6  7  8  9
%                   3 10 11  6 12 13  9 14 15 
%                   9 14 15  8 16 17  7 18 19 ];

Xglobal = zeros(3,(N+1)^2);
Yglobal = zeros(3,(N+1)^2);
dXdXiglobal  = zeros(3,(N+1)^2);
dXdEtaglobal = zeros(3,(N+1)^2);
dYdXiglobal  = zeros(3,(N+1)^2);
dYdEtaglobal = zeros(3,(N+1)^2);
for i=1:3

[Xglobal(i,:),dXdXiglobal(i,:),dXdEtaglobal(i,:)] = ...
                              transfinitemapping(xglobal(corners(i,:)));
[Yglobal(i,:),dYdXiglobal(i,:),dYdEtaglobal(i,:)] = ...
                              transfinitemapping(yglobal(corners(i,:)));

pcolor(reshape(Xglobal(i,:),N+1,N+1),reshape(Yglobal(i,:),N+1,N+1),(i-1)/2*ones(N+1))
hold on
end

Jglobal = dXdXiglobal.*dYdEtaglobal-dXdEtaglobal.*dYdXiglobal;


plot(xglobal,yglobal,'o','markerface','b')
axis equal
axis([-1.2 1.2 -.6 1.2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%