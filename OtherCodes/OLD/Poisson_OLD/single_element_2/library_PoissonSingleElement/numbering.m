function [] = numbering()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% [globalnr_0,globalnr_1h,globalnr_1v,globalnr_2]                         %
%                                      = numbering(N,numColumns,numRows)  %
%                                                                         %
% Global numbering for single grid problems for single element in 2D      %
%                                                                         %
% Global numbering for point values (0-cells)                             %
% Global numbering for lines values (1-cells)                             %
% Global numbering for surface values (2-cells)                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_0(1:N+1,1:N+1) = reshape(1:(N+1)^2,N+1,N+1);

% disp('  '); disp('globalnr_0 = '); disp('  ');
% disp(num2str(flipud(globalnr_0')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1v(1:N+1,1:N) = reshape(1:N*(N+1),N+1,N);
globalnr_1h(1:N,1:N+1) = N*(N+1) + reshape(1:N*(N+1),N+1,N)';

% disp('  '); disp('globalnr_1v = '); disp('  ');
% disp(num2str(flipud(globalnr_1v')));
% disp('  '); disp('globalnr_1h = '); disp('  ');
% disp(num2str(flipud(globalnr_1h')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2(1:N,1:N) = reshape(1:N^2,N,N);

% disp('  '); disp('globalnr_2 = '); disp('  ');
% disp(num2str(flipud(globalnr_2')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = globalnr_0(end,end);
nr_1 = globalnr_1h(end,end);
nr_2 = globalnr_2(end,end);