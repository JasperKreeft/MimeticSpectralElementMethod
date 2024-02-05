function numbering_adaptive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% [] = numbering(varargin)                                                %
%                                                                         %
% Global numbering for single grid problems in 2D                         %
%                                                                         %
% Global numbering for point values (0-cells)                             %
% Global numbering for lines values (1-cells)                             %
% Global numbering for surface values (2-cells)                           %
%                                                                         %
% numbering('single') is single element numbering                         %
% numbering('multi')  is multi  element numbering                         %
% numbering()         is multi  element numbering                         %
%                                                                         %
% global output: globalnr_0 globalnr_1v globalnr_1h globalnr_2            %
%                nr_0 nr_1 nr_2                                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 27-06-2011                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Nadaptive numElements
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

Nmax = max(Nadaptive);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_0 = zeros((Nmax+1)^2,numElements);

ind0 = 0;
for i=1:numElements
    ind1 = 1:(Nadaptive(i)+1)^2;
    globalnr_0(ind1,i) = ind0 + ind1;
    ind0 = globalnr_0(ind1(end),i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1h = zeros(Nmax*(Nmax+1),numElements);
globalnr_1v = zeros(Nmax*(Nmax+1),numElements);

ind0 = 0;
for i=1:numElements
    ind1 = 1:Nadaptive(i)*(Nadaptive(i)+1);
    globalnr_1v(ind1,i) = ind0 + ind1;
    ind0 = globalnr_1v(ind1(end),i);
    globalnr_1h(ind1,i) = ind0 + ind1;
    ind0 = globalnr_1h(ind1(end),i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2 = zeros(Nmax^2,numElements);
ind = 0;
for i=1:numElements
    ind2 = 1:Nadaptive(i)^2;
    globalnr_2(ind2,i) = ind+(1:Nadaptive(i)^2);
    ind = globalnr_2(Nadaptive(i)^2,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = max(globalnr_0(:,end));
nr_1 = max(globalnr_1h(:,end));
nr_2 = max(globalnr_2(:,end));