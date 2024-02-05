function [A b nodesOut firstVal dirichletValues] = bcZeroFormTwoD(BCs,A,b,boundaryNodes,boundaryValues)
%
%   bcZeroFormTwoD applies the boundary conditions for zero forms in 2D
%
%   [A b nodesOut firstVal dirichletValues] =
%       bcZeroFormTwoD(BCs,A,b,boundaryNodes,boundaryValues)
%
%   input:
%       BCs             :: boundary conditions
%       A               :: system matrix
%       B               :: system rhs
%       boundaryNodes   :: boundary nodes numbering
%       boundaryValues  :: boundary nodes values
% 
%   output:
%       A               :: updated system matrix
%       b               :: updated system rhs
%       nodesOut        :: global numbering of dirichlet boundary nodes only
%       firstVal        :: index of first occurences of each unique dirichlet boundary node 
%       dirichletValues :: values at dirichlet boundary nodes
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 26/10/2011 $  

%-------------------------------------------------------------------------%
% check bcs                                                               %
%-------------------------------------------------------------------------%

    dirichletSides = find(BCs); % sides with dirichlet bc
    neumannSides = find(BCs-1); % sides with neumann bc
    dirichletFlag = 1;          % dirichlet flag
    if isempty(dirichletSides)  % if empty then...
        dirichletFlag = 0;      % ...no dirichlet bc
    end
    neumannFlag = 1;            % dirichletflag
    if isempty(neumannSides)    % if empty then...
        neumannFlag = 0;        % ...no neumann bc
    end

%-------------------------------------------------------------------------%
% apply neumann                                                           %
%-------------------------------------------------------------------------%   
    
    if neumannFlag
           
%         for nN = neumannSides
%             neumannNodes = boundaryNodes{nN};                   % numbering
%             neumannValues = boundaryValues{nN};                 % values
%             b(neumannNodes) = b(neumannNodes)+neumannValues;	% rhs contribution
%         end
%          
%         %%% output
%         nodesOut = [];          % not used in case of neumann
%         firstVal = [];          % not used in case of neumann
%         dirichletValues = [];   % not used in case of neumann
        
    end
 
%-------------------------------------------------------------------------%
% apply dirichlet                                                         %
%-------------------------------------------------------------------------%       
    
    if dirichletFlag
        
        %%% numbering and values
        dirichletNodes = [];
        dirichletValues = [];
        for nD = 1:length(dirichletSides) 
            dNodes = reshape(boundaryNodes{nD},numel(boundaryNodes{nD}),1);
            dirichletNodes = vertcat(dirichletNodes,dNodes);
            dValues = reshape(boundaryValues{nD},numel(boundaryValues{nD}),1);
            dirichletValues = vertcat(dirichletValues,dValues);            
        end
        
        %%% sort and remove duplicates
        [dirichletNodes firstVal newPos] = unique(dirichletNodes,'first'); %#ok<NASGU>
        % note:
        % dirichletNodes = gn of dirichlet nodes (without duplicates)
        % firstVal = index of first occurence of each unique value 
        % newPos = updated position in newly ordered arrangement
              
        %%% update system rhs
        b = b-A(:,dirichletNodes)*dirichletValues(firstVal); 

        %%% update system matrix
        A(dirichletNodes,:) = [];                   
        A(:,dirichletNodes) = [];                    
        b(dirichletNodes) = [];     

        %%% output
        nodesOut = dirichletNodes;
                
    end
        
end
    