function gn = GlobalNumberingMomentumPrimal(n,p,varargin)
%GlobalNumberingTwoFormPrimal computes the global numbering of 2-form
%degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingTwoFormPrimal(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the momentum finite-volumes of the mesh. gn is a matrix
%   N x p^2, where N is the total number of elements.
%
%   First all U/V (depending on orientation of primal mesh) finite-volumes are numbered, and then all V/U finite-volumes.
%
%
%   See also: GLOBALNUMBERINGZEROFORMDUAL, PLOTGLOBALNUMBERINZEROFORMPRIMAL

%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2012/03/02 $
    
    if size(varargin,2)
        periodic = varargin{1};
    else
        periodic = [false false];
    end

    if length(n)>2
        % check if n has more than two values
        disp(sprintf(':: n :: has too many values'));
        return
    elseif length(n)==1
            n = [n n];
    end
    
    if prod(n)==0
        % check if n has a zero value
        disp(sprintf(':: n :: has at least one zero'));
        return
    end
    
    % compute the number of elements
    nElements = n(1)*n(2);

    % compute the total number of degrees of freedom
    dof = p*(p+1)*nElements;
    
    % compute the global numbering
%     gn = [reshape(uint32(1):uint32(dof), p*(p+1), nElements)'   reshape(uint32(1+dof):uint32(2*dof), p*(p+1), nElements)'];
    gn.Xi = reshape(uint32(1):uint32(dof), p*(p+1), nElements)';
    gn.Eta = gn.Xi;
    
    %% Global numbering of momenta
    
    % This yields the exact same number of finite-volumes as velocity
    % 1-forms.
    
    % Global numbering of 1-forms
    gnOne = GlobalNumberingOneFormPrimalPeriodic(n,p,periodic);
    dofOne = 0.5*size(gnOne,2);
    gnOneXi = gnOne(:,1:dofOne);
    gnOneEta = gnOne(:,(dofOne+1):end);
    
    %%% XI
    count = 1;
    temp1 = gnOneXi';
    temp1 = double(temp1(:));
    temp2 = zeros(size(temp1,1),2);
    for resultNum = 1:size(temp1,1)
        
        [new,indexRepeat] = min(abs(temp2(:,2)-temp1(resultNum)));
        if (new)
            temp2(resultNum,:) = [count temp1(resultNum)];
            count = count+1;
        else
            temp2(resultNum,:) = [temp2(indexRepeat,1) temp1(resultNum)];
        end
    end
    gn.XiG = (reshape(temp2(:,1),size(temp2,1)/nElements,nElements))';
    clear temp1 temp2 resultNum count indexRepeat;
    
    %%% ETA
    count = 1;
    temp1 = gnOneEta';
    temp1 = double(temp1(:));
    temp2 = zeros(size(temp1,1),2);
    for resultNum = 1:size(temp1,1)
        
        [new,indexRepeat] = min(abs(temp2(:,2)-temp1(resultNum)));
        if (new)
            temp2(resultNum,:) = [count temp1(resultNum)];
            count = count+1;
        else
            temp2(resultNum,:) = [temp2(indexRepeat,1) temp1(resultNum)];
        end
    end
    gn.EtaG = (reshape(temp2(:,1),size(temp2,1)/nElements,nElements))';
    

end