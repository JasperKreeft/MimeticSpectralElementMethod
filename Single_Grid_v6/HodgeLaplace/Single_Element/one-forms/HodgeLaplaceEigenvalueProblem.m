%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Single-element Hodge-Laplace on 1-forms problem                         %
% Hodge Laplace is curl-curl + div-div                                    %
%                                                                         %
% written by Jasper Kreeft (2010)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N w
global h e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Problem = "DivDiv"; %"CurlCurl"; % "HodgeLaplacian"; %  

NrCellRange = 2:2:20;

cc = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

error = zeros(10); er = 0;

for N=NrCellRange
    
disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

[xi,w] = GLLnodes(N);     eta = xi;                % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square('CurlCurl',cc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgrad(N);

D = div(N);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system

switch Problem
    case "HodgeLaplacian"
        L = M1*NG/M0*NG'*M1+D'*M2*D;
        % L = -[ M0 NG'*M1
        %       M1*NG -D'*M2*D ];
        case "CurlCurl"
        L = M1*NG/M0*NG'*M1;
        % L = -[ M0 NG'*M1
        %       M1*NG eps*ones(2*N*(N+1)) ];
    case "DivDiv"
        L = D'*M2*D;
end

E = sort(eig(full(L),full(M1)));

% RHS = [ zeros((N+1)^2) zeros((N+1)^2,2*N*(N+1))
%         zeros(2*N*(N+1),(N+1)^2)  M1            ];
% 
% E = sort(eig(full(L),full(RHS)));

E(abs(E)<.2)=[];

switch Problem
    case "HodgeLaplacian"
        exact = [1 1 2 2 4 4 5 5 5 5]'; % 8 8 9 9 10 10 10 10 13 13
    case "CurlCurl"
        exact = [1 1 2 4 4 5 5 8 9 9]';
    case "DivDiv"
        exact = [2 5 5 8 10 10 13 13 17 17]';
end

nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));
disp([0;E(1:min(length(E),20))])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

plot_eigenvalues(E)

plot_convergence_eigenvalues("N",NrCellRange,error,exact)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%