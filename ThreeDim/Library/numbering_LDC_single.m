function []=numbering_LDC_single()

global N
global globalnr_0
global globalnr_1x globalnr_1y globalnr_1z
global globalnr_2x globalnr_2y globalnr_2z
global globalnr_3
global nr_0 nr_1 nr_2 nr_3

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_0 = zeros((N+1)^3,1);

el1 = reshape((1:(N+1)^3),N+1,N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^3,1);
nr_0 = max(globalnr_0(:,1));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1x = zeros(N*(N+1)^2,1);
globalnr_1y = zeros(N*(N+1)^2,1);
globalnr_1z = zeros(N*(N+1)^2,1);

el1z = reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1z(:,1) = reshape(el1z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,1));

el1y = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1y(:,1) = reshape(el1y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,1));

el1x = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1x(:,1) = reshape(el1x,N*(N+1)^2,1);
nr_1 = max(globalnr_1x(:,1));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2x = zeros(N^2*(N+1),1);
globalnr_2y = zeros(N^2*(N+1),1);
globalnr_2z = zeros(N^2*(N+1),1);

el1x = reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2x(:,1) = reshape(el1x,N^2*(N+1),1);
ind = max(globalnr_2x(:,1));

el1y = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2y(:,1) = reshape(el1y,N^2*(N+1),1);
ind = max(globalnr_2y(:,1));

el1z = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2z(:,1) = reshape(el1z,N^2*(N+1),1);
nr_2 = max(globalnr_2z(:,1));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_3 = zeros(N^3,1);

globalnr_3(:,1) = (1:N^3);
nr_3 = globalnr_3(N^3,1);
