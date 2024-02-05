function [E10 E21] = incidenceMatrix(nElement,p)
%
%   incidenceMatrix produces the incidence matrix for a domain consisting
%   of nElement elements of order p (it combines xi en eta contributions)
%
%   incidenceMatrix(nElement,p)
%
%   input:
%       nElement	:: number of elements [nElement_x nElemen_y]
%       p           :: element order [p_x p_y]
%
%   output:
%       E10         :: incidence matrix relating one forms to zero forms
%       E21         :: incidence matrix relating two forms to one forms 
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 01/12/2011 $

%-------------------------------------------------------------------------%
% input checks                                                            %
%-------------------------------------------------------------------------%   
    
    % check nElement
    if length(nElement)>2
        fprintf(':: nElement :: has too many values \n');
        return
    elseif length(nElement)==1
        nElement = [nElement nElement];
    end    
    if prod(nElement)==0
        fprintf(':: nElement :: has at least one zero \n');
        return
    end 
    % check p
    if length(p)>2
        fprintf(':: p :: has too many values \n');
        return
    elseif length(p)==1
        p = [p p];
    end    
    if prod(p)==0
        fprintf(':: p :: has at least one zero \n');
        return
    end 

%-------------------------------------------------------------------------%
% generate k-forms (k=0,1,2)                                              %
%-------------------------------------------------------------------------%   

    zeroForms = GlobalNumberingZeroFormTwoD(nElement,p);
    oneForms = GlobalNumberingOneFormTwoD(nElement,p);
    twoForms = GlobalNumberingTwoFormTwoD(nElement,p);

%-------------------------------------------------------------------------%
% E10                                                                     %
%-------------------------------------------------------------------------%  

    E10 = zeros(oneForms(end,end),max(max(zeroForms)));
    E10x = full(dZeroXiTwoD(p));
    E10y = full(dZeroEtaTwoD(p));
    nOneFormsX = p(1)*(p(2)+1);
    % nOneFormsY = p(2)*(p(1)+1);
    
    for element = 1:prod(nElement)
        E10(oneForms(element,1:nOneFormsX),zeroForms(element,:))= E10x; 
        E10(oneForms(element,(nOneFormsX+1):end),zeroForms(element,:))= E10y;  
    end      

%-------------------------------------------------------------------------%
% E10                                                                     %
%-------------------------------------------------------------------------%      
    
    E21 = zeros(oneForms(end,end),twoForms(end,end));

    E21x = full(dOneXiTwoD(p));
    E21y = full(dOneEtaTwoD(p));
    for element=1:prod(nElement)
        E21(oneForms(element,1:nOneFormsX),twoForms(element,:))= E21x';
        E21(oneForms(element,(nOneFormsX+1):end),twoForms(element,:))= E21y';    
    end
    E21 = E21';
    
end
