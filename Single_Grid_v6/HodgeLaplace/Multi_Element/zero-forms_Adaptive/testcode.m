clear all
close all
clc

global Nadaptive numElements Nmax
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

Nadaptive = [ 2 4 ; 1 2 ];

numElements = numel(Nadaptive);
numRows    = size(Nadaptive,1);
numColumns = size(Nadaptive,2);

Nmax = max(max(Nadaptive));

numbering_adaptive

Na = zeros(numColumns*(numRows-1)+numRows*(numColumns-1),2);
for c=1:numColumns-1
    rc = (1:numRows)+(c-1)*numRows;
    Na(rc,:) = Nadaptive(:,c:c+1);
end
for r=1:numRows-1
    rc = (1:numColumns) + (r-1)*numColumns + numRows*(numColumns-1);
    Na(rc,:) = Nadaptive(r:r+1,:)';
end
% Na = unique(Na,'rows');


globalnr0 = 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Connectivity

nr_0_in_element(i) = (Nadaptive+1).^2;
for i=1:numElements
    for b=1:4
        dif_N = Nadaptive(i) - Nadaptive(cnectN(b));
        nr_0_in_element(i) = nr_0_in_element(i) - (dif_N>0)*dif_N;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




nr_0_in_element = zeros(numElements,1);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
    
        if rc==1 && numColumns>1
            nr_0_in_element(rc) = Nadaptive(r,c)^2 + ...
                                 (Nadaptive(r,c)<=Nadaptive(r,c+1))*Nadaptive(r,c) + ...
                                 (Nadaptive(r,c)<=Nadaptive(r+1,c))*Nadaptive(r,c);

        elseif c<numColumns


        elseif c==numColumns && r==1
            nr_0_in_element(rc) = (Nadaptive(r,c)-1)^2 + ...
                                 (Nadaptive(r,c)<=Nadaptive(r,c-1))*Nadaptive(r,c) + ...
                                 (Nadaptive(r,c)<=Nadaptive(r+1,c))*Nadaptive(r,c);
        end

    end
end




nr_0_in_element




