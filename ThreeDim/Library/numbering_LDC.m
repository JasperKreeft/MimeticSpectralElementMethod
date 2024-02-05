function []=numbering_LDC()

global N
global globalnr_0
global globalnr_1x globalnr_1y globalnr_1z
global globalnr_2x globalnr_2y globalnr_2z
global globalnr_3
global nr_0 nr_1 nr_2 nr_3

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

el0 = zeros(N+1,N+1,N+1);

globalnr_0 = zeros((N+1)^3,8);

el1 = reshape((1:(N+1)^3),N+1,N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^3,1);
ind = max(globalnr_0(:,1));

el2 = el0;
el2(1,:,:) = el1(N+1,:,:);
el2(2:N+1,:,:) = reshape(ind+(1:N*(N+1)*(N+1)),N,N+1,N+1);
globalnr_0(:,2) = reshape(el2,(N+1)^3,1);
ind = max(globalnr_0(:,2));

el3 = el0;
el3(:,1,:) = el1(:,N+1,:);
el3(:,2:N+1,:) = reshape(ind+(1:N*(N+1)*(N+1)),N+1,N,N+1);
globalnr_0(:,3) = reshape(el3,(N+1)^3,1);
ind = max(globalnr_0(:,3));

el4 = el0;
el4(:,1,:) = el2(:,N+1,:);
el4(1,:,:) = el3(N+1,:,:);
el4(2:N+1,2:N+1,:) = reshape(ind+(1:N*N*(N+1)),N,N,N+1);
globalnr_0(:,4) = reshape(el4,(N+1)^3,1);
ind = max(globalnr_0(:,4));

el5 = el0;
el5(:,:,1) = el1(:,:,N+1);
el5(:,:,2:N+1) = reshape(ind+(1:N*(N+1)*(N+1)),N+1,N+1,N);
globalnr_0(:,5) = reshape(el5,(N+1)^3,1);
ind = max(globalnr_0(:,5));

el6 = el0;
el6(:,:,1) = el2(:,:,N+1);
el6(1,:,:) = el5(N+1,:,:);
el6(2:N+1,:,2:N+1) = reshape(ind+(1:N*N*(N+1)),N,N+1,N);
globalnr_0(:,6) = reshape(el6,(N+1)^3,1);
ind = max(globalnr_0(:,6));

el7 = el0;
el7(:,:,1) = el3(:,:,N+1);
el7(:,1,:) = el5(:,N+1,:);
el7(:,2:N+1,2:N+1) = reshape(ind+(1:N*N*(N+1)),N+1,N,N);
globalnr_0(:,7) = reshape(el7,(N+1)^3,1);
ind = max(globalnr_0(:,7));

% 8th element
el8 = el0;
el8(:,:,1) = el4(:,:,N+1);
el8(:,1,:) = el6(:,N+1,:);
el8(1,:,:) = el7(N+1,:,:);
el8(2:N+1,2:N+1,2:N+1) = reshape(ind+(1:N^3),N,N,N);
globalnr_0(:,8) = reshape(el8,(N+1)^3,1);
nr_0 = max(globalnr_0(:,8));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1x = zeros(N*(N+1)^2,8);
globalnr_1y = zeros(N*(N+1)^2,8);
globalnr_1z = zeros(N*(N+1)^2,8);

el0x = zeros(N,N+1,N+1);
el0y = zeros(N,N+1,N+1);
el0z = zeros(N,N+1,N+1);

el1z = reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1z(:,1) = reshape(el1z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,1));

el1y = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1y(:,1) = reshape(el1y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,1));

el1x = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1x(:,1) = reshape(el1x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,1));

el2z = el0z;
el2z(:,1,:) = el1z(:,N+1,:);
el2z(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1z(:,2) = reshape(el2z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,2));

el2y = el0y;
el2y(:,1,:) = el1y(:,N+1,:);
el2y(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1y(:,2) = reshape(el2y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,2));

el2x = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1x(:,2) = reshape(el2x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,2));

el3z = el0z;
el3z(:,:,1) = el1z(:,:,N+1);
el3z(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1z(:,3) = reshape(el3z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,3));

el3y = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1y(:,3) = reshape(el3y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,3));

el3x = el0x;
el3x(:,1,:) = el1x(:,N+1,:);
el3x(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1x(:,3) = reshape(el3x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,3));

el4z = el0z;
el4z(:,:,1) = el2z(:,:,N+1);
el4z(:,1,:) = el3z(:,N+1,:);
el4z(:,2:N+1,2:N+1) = ind+reshape(1:N^3,N,N,N);
globalnr_1z(:,4) = reshape(el4z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,4));

el4y = el0y;
el4y(:,1,:) = el3y(:,N+1,:);
el4y(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1y(:,4) = reshape(el4y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,4));

el4x = el0x;
el4x(:,1,:) = el2x(:,N+1,:);
el4x(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1x(:,4) = reshape(el4x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,4));

el5z = ind+reshape(1:N*(N+1)^2,N,N+1,N+1);
globalnr_1z(:,5) = reshape(el5z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,5));

el5y = el0y;
el5y(:,:,1) = el1y(:,:,N+1);
el5y(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1y(:,5) = reshape(el5y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,5));

el5x = el0x;
el5x(:,:,1) = el1x(:,:,N+1);
el5x(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1x(:,5) = reshape(el5x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,5));

el6z = el0z;
el6z(:,1,:) = el5z(:,N+1,:);
el6z(:,2:N+1,:) = ind+reshape(1:N^2*(N+1),N,N,N+1);
globalnr_1z(:,6) = reshape(el6z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,6));

el6y = el0y;
el6y(:,:,1) = el2y(:,:,N+1);
el6y(:,1,:) = el5y(:,N+1,:);
el6y(:,2:N+1,2:N+1) = ind+reshape(1:N^3,N,N,N);
globalnr_1y(:,6) = reshape(el6y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,6));

el6x = el0x;
el6x(:,:,1) = el2x(:,:,N+1);
el6x(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1x(:,6) = reshape(el6x,N*(N+1)^2,1);
ind = max(globalnr_1x(:,6));

el7z = el0z;
el7z(:,:,1) = el5z(:,:,N+1);
el7z(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1z(:,7) = reshape(el7z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,7));

el7y = el0y;
el7y(:,:,1) = el3y(:,:,N+1);
el7y(:,:,2:N+1) = ind+reshape(1:N^2*(N+1),N,N+1,N);
globalnr_1y(:,7) = reshape(el7y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,7));

el7x = el0x;
el7x(:,:,1) = el3x(:,:,N+1);
el7x(:,1,:) = el5x(:,N+1,:);
el7x(:,2:N+1,2:N+1) = ind+reshape(1:N^3,N,N,N);
globalnr_1x(:,7) = reshape(el7x,N*(N+1)^2,1);
ind= max(globalnr_1x(:,7));

% 8th element
el8z = el0z;
el8z(:,:,1) = el6z(:,:,N+1);
el8z(:,1,:) = el7z(:,N+1,:);
el8z(:,2:N+1,2:N+1) = reshape(ind+(1:N^3),N,N,N);
globalnr_1z(:,8) = reshape(el8z,N*(N+1)^2,1);
ind = max(globalnr_1z(:,8));

el8y = el0y;
el8y(:,:,1) = el4y(:,:,N+1);
el8y(:,1,:) = el7y(:,N+1,:);
el8y(:,2:N+1,2:N+1) = reshape(ind+(1:N^3),N,N,N);
globalnr_1y(:,8) = reshape(el8y,N*(N+1)^2,1);
ind = max(globalnr_1y(:,8));

el8x = el0x;
el8x(:,:,1) = el4x(:,:,N+1);
el8x(:,1,:) = el6x(:,N+1,:);
el8x(:,2:N+1,2:N+1) = reshape(ind+(1:N^3),N,N,N);
globalnr_1x(:,8) = reshape(el8x,N*(N+1)^2,1);
nr_1= max(globalnr_1x(:,8));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2x = zeros(N^2*(N+1),8);
globalnr_2y = zeros(N^2*(N+1),8);
globalnr_2z = zeros(N^2*(N+1),8);

el0x = zeros(N+1,N,N);
el0y = zeros(N+1,N,N);
el0z = zeros(N+1,N,N);

el1x = reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2x(:,1) = reshape(el1x,N^2*(N+1),1);
ind = max(globalnr_2x(:,1));

el1y = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2y(:,1) = reshape(el1y,N^2*(N+1),1);
ind = max(globalnr_2y(:,1));

el1z = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2z(:,1) = reshape(el1z,N^2*(N+1),1);
ind = max(globalnr_2z(:,1));

el2x = el0x;
el2x(1,:,:) = el1x(N+1,:,:);
el2x(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2x(:,2) = reshape(el2x,N^2*(N+1),1);
ind = max(globalnr_2x(:,2));

el2y = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2y(:,2) = reshape(el2y,N^2*(N+1),1);
ind = max(globalnr_2y(:,2));

el2z = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2z(:,2) = reshape(el2z,N^2*(N+1),1);
ind = max(globalnr_2z(:,2));

el3x = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2x(:,3) = reshape(el3x,N^2*(N+1),1);
ind = max(globalnr_2x(:,3));

el3y = el0y;
el3y(1,:,:) = el1y(N+1,:,:);
el3y(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2y(:,3) = reshape(el3y,N^2*(N+1),1);
ind = max(globalnr_2y(:,3));

el3z = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2z(:,3) = reshape(el3z,N^2*(N+1),1);
ind = max(globalnr_2z(:,3));

el4x = el0x;
el4x(1,:,:) = el3x(N+1,:,:);
el4x(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2x(:,4) = reshape(el4x,N^2*(N+1),1);
ind = max(globalnr_2x(:,4));

el4y = el0y;
el4y(1,:,:) = el2y(N+1,:,:);
el4y(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2y(:,4) = reshape(el4y,N^2*(N+1),1);
ind = max(globalnr_2y(:,4));

el4z = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2z(:,4) = reshape(el4z,N^2*(N+1),1);
ind = max(globalnr_2z(:,4));

el5x = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2x(:,5) = reshape(el5x,N^2*(N+1),1);
ind = max(globalnr_2x(:,5));

el5y = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2y(:,5) = reshape(el5y,N^2*(N+1),1);
ind = max(globalnr_2y(:,5));

el5z = el0z;
el5z(1,:,:) = el1z(N+1,:,:);
el5z(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2z(:,5) = reshape(el5z,N^2*(N+1),1);
ind = max(globalnr_2z(:,5));

el6x = el0x;
el6x(1,:,:) = el5x(N+1,:,:);
el6x(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2x(:,6) = reshape(el6x,N^2*(N+1),1);
ind = max(globalnr_2x(:,6));

el6y = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2y(:,6) = reshape(el6y,N^2*(N+1),1);
ind = max(globalnr_2y(:,6));

el6z = el0z;
el6z(1,:,:) = el2z(N+1,:,:);
el6z(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2z(:,6) = reshape(el6z,N^2*(N+1),1);
ind = max(globalnr_2z(:,6));

el7x = ind+reshape(1:N^2*(N+1),N+1,N,N);
globalnr_2x(:,7) = reshape(el7x,N^2*(N+1),1);
ind = max(globalnr_2x(:,7));

el7y = el0y;
el7y(1,:,:) = el5y(N+1,:,:);
el7y(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2y(:,7) = reshape(el7y,N^2*(N+1),1);
ind = max(globalnr_2y(:,7));

el7z = el0z;
el7z(1,:,:) = el3z(N+1,:,:);
el7z(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2z(:,7) = reshape(el7z,N^2*(N+1),1);
ind = max(globalnr_2z(:,7));


% 8th element
el8x = el0x;
el8x(1,:,:) = el7x(N+1,:,:);
el8x(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2x(:,8) = reshape(el8x,N^2*(N+1),1);
ind = max(globalnr_2x(:,8));

el8y = el0y;
el8y(1,:,:) = el6y(N+1,:,:);
el8y(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2y(:,8) = reshape(el8y,N^2*(N+1),1);
ind = max(globalnr_2y(:,8));

el8z = el0z;
el8z(1,:,:) = el4z(N+1,:,:);
el8z(2:N+1,:,:) = ind+reshape(1:N^3,N,N,N);
globalnr_2z(:,8) = reshape(el8z,N^2*(N+1),1);
nr_2 = max(globalnr_2z(:,8));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_3 = zeros(N^3,8);
for i=1:8
    globalnr_3(:,i) = (i-1)*N^3+(1:N^3);
end
nr_3 = globalnr_3(N^3,8);
