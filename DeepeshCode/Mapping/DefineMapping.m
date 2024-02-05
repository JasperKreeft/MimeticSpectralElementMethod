function [phi g g11 g12 g22 dPhiXdXi dPhiYdXi dPhiXdEta dPhiYdEta varargout] = DefineMapping(n,mapXi_Coeff1,mapXi_Coeff2,mapEta_Coeff1,mapEta_Coeff2,deltaX,deltaY,map,curved,varargin)

    %% Memory allocation
    phi = cell(n(1)*n(2),1);% allocate memory space for the mapping

    g11 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g12 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g22 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    g = cell(n(1)*n(2),1);% allocate memory space for the mapping
    
    dPhiXdXi = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiXdEta = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiYdXi = cell(n(1)*n(2),1);% allocate memory space for the mapping
    dPhiYdEta = cell(n(1)*n(2),1);% allocate memory space for the mapping
    
    if (size(varargin,2))        
        covariantMetric = varargin{1};
        vectorTransformation = varargin{2};
    end
    
    if (covariantMetric)
        gU11 = cell(n(1)*n(2),1);% allocate memory space for the mapping
        gU12 = cell(n(1)*n(2),1);% allocate memory space for the mapping
        gU22 = cell(n(1)*n(2),1);% allocate memory space for the mapping
    end
    
    if (vectorTransformation)
        dPhiXidX = cell(n(1)*n(2),1);% allocate memory space for the mapping
        dPhiEtadX = cell(n(1)*n(2),1);% allocate memory space for the mapping
        dPhiXidY = cell(n(1)*n(2),1);% allocate memory space for the mapping
        dPhiEtadY = cell(n(1)*n(2),1);% allocate memory space for the mapping
    end
    


    %% Assign mappings
    
    if (strcmp(map,'Normal')) % Normal (scaling in xi and eta directions) map
        % If 'curved' variable is missing, set it to zero
        if ~exist('curved','var')
            curved = 0;
        end
        
        % loop over the elements and generate the mappings
        for elementXi = 1:n(1)
            for elementEta = 1:n(2)
                element = (elementXi-1)*n(2) + elementEta;
                phi{element} = @(xi,eta) (deal(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi)+curved*sin(pi*xi).*sin(pi*eta), ... 
                                                mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)+curved*sin(pi*xi).*sin(pi*eta)));
                % Jacobian elements
                dPhiXdXi{element} = @(xi,eta) mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta);
                dPhiXdEta{element} = @(xi,eta) zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta);
                dPhiYdXi{element} = @(xi,eta) zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta);
                dPhiYdEta{element} = @(xi,eta) mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta);

                % sqrt of determinant of metric tensor
    %             g{element} = @(xi,eta) abs((mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).*(mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) ...
    %                  - (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).*(zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta)));
                g{element} = @(xi,eta) abs(dPhiXdXi{element}(xi,eta).*dPhiYdEta{element}(xi,eta) - dPhiYdXi{element}(xi,eta).*dPhiXdEta{element}(xi,eta));
                g11{element} = @(xi,eta) (dPhiXdEta{element}(xi,eta).^2 + dPhiYdEta{element}(xi,eta).^2)./(g{element}(xi,eta).^2);
                g12{element} = @(xi,eta) (-dPhiXdEta{element}(xi,eta).*dPhiXdXi{element}(xi,eta) -dPhiYdEta{element}(xi,eta).*dPhiYdXi{element}(xi,eta))./(g{element}(xi,eta).^2);
                g22{element} = @(xi,eta) (dPhiXdXi{element}(xi,eta).^2 + dPhiYdXi{element}(xi,eta).^2)./(g{element}(xi,eta).^2);
            end
        end

        for element = 1:n(1)*n(2)
            
            % g^{11}
%             g11{element} = @(xi,eta) ((zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).^2 + ...
%                                 (mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).^2)./ ... 
%                                 (((mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).*(mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) ...
%                  - (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).*(zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta))).^2);
%             g11{element} = @(xi,eta) (dPhiXdEta{element}(xi,eta).^2 + dPhiYdEta{element}(xi,eta).^2)./(g{element}(xi,eta).^2);
            % g^{12} = g^{21}
%             g12{element} = @(xi,eta) ((mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).* ... 
%                                         (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) + ...
%                                             (zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta)).* ... 
%                                                 (mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)))./ ...
%                                                 (((mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).*(mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) ...
%                  - (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).*(zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta))).^2);
%             g12{element} = @(xi,eta) (-dPhiXdEta{element}(xi,eta).*dPhiXdXi{element}(xi,eta) -dPhiYdEta{element}(xi,eta).*dPhiYdXi{element}(xi,eta))./(g{element}(xi,eta).^2);
            % g^{22}
%             g22{element} = @(xi,eta) ((mapXi_Coeff1(elementXi)*ones(size(xi)) + curved*pi*cos(pi*xi).*sin(pi*eta)).^2+(zeros(size(xi)) ...
%                             + curved*pi*cos(pi*xi).*sin(pi*eta)).^2)./(((mapXi_Coeff1(elementXi)*ones(size(xi)) + ...
%                             curved*pi*cos(pi*xi).*sin(pi*eta)).*(mapEta_Coeff1(elementEta)*ones(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)) ...
%                  - (zeros(size(xi)) + curved*pi*sin(pi*xi).*cos(pi*eta)).*(zeros(size(xi))  + curved*pi*cos(pi*xi).*sin(pi*eta))).^2);
%             g22{element} = @(xi,eta) (dPhiXdXi{element}(xi,eta).^2 + dPhiYdXi{element}(xi,eta).^2)./(g{element}(xi,eta).^2);
                
            if (covariantMetric)
                gU11{element} = @(xi,eta) (dPhiXdXi{element}(xi,eta).^2 + dPhiYdXi{element}(xi,eta).^2);
                gU12{element} = @(xi,eta) (dPhiXdEta{element}(xi,eta).*dPhiXdXi{element}(xi,eta)+dPhiYdEta{element}(xi,eta).*dPhiYdXi{element}(xi,eta));
                gU22{element} = @(xi,eta) (dPhiXdEta{element}(xi,eta).^2 + dPhiYdEta{element}(xi,eta).^2);
                
                if (vectorTransformation)
                    
                    % (mapXi_Coeff1(elementXi)*ones(size(xi)) + pi*curved*cos(pi*xi).*sin(pi*eta)).*(mapEta_Coeff1(elementEta)*ones(size(eta)) + pi*curved*sin(pi*xi).*cos(pi*eta)) ...
                    % - (pi*curved*cos(pi*xi).*sin(pi*eta)).*(pi*curved*sin(pi*xi).*cos(pi*eta))
                    
                    % Jacobian elements
                    dPhiXidX{element} = @(xi,eta) (dPhiYdEta{element}(xi,eta))./ (g{element}(xi,eta));
                    dPhiEtadX{element} = @(xi,eta) (-dPhiYdXi{element}(xi,eta))./ (g{element}(xi,eta));
                    dPhiXidY{element} = @(xi,eta) (-dPhiXdEta{element}(xi,eta))./ (g{element}(xi,eta));
                    dPhiEtadY{element} = @(xi,eta) (dPhiXdXi{element}(xi,eta))./ (g{element}(xi,eta));
                end

            end             
            
           
        end
        
        if (covariantMetric)
            varargout{1} = gU11;
            varargout{2} = gU12;
            varargout{3} = gU22;
        end
        
        if(vectorTransformation)
            varargout{4} = dPhiXidX;
            varargout{5} = dPhiEtadX;
            varargout{6} = dPhiXidY;
            varargout{7} = dPhiEtadY;
        end
        
        
    elseif (strcmp(map,'HalfDisc'))
        
        % Mapping
        for elementXi = 1:n(1)
            for elementEta = 1:n(2)
                element = (elementXi-1)*n(2) + elementEta;
                phi{element} = @(xi,eta) deal((mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi)), ...
                    (mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi)));
            end
        end
    
        % Jacobian elements and Metric tensor components
        for elementXi = 1:n(1)
            for elementEta = 1:n(2)

                element = (elementXi-1)*n(2) + elementEta;

                % Jacobian elements
                dPhiXdXi{element} = @(xi,eta) -mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi));
                dPhiXdEta{element} = @(xi,eta) mapEta_Coeff1(elementEta)*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi));
                dPhiYdXi{element} = @(xi,eta) mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi));
                dPhiYdEta{element} = @(xi,eta) mapEta_Coeff1(elementEta)*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi));
                
                % Metric Tensor components etc.
                g{element} = @(xi,eta) abs(dPhiXdXi{element}(xi,eta).*dPhiYdEta{element}(xi,eta) - dPhiYdXi{element}(xi,eta).*dPhiXdEta{element}(xi,eta));
                g11{element} = @(xi,eta) (dPhiXdEta{element}(xi,eta).^2 + dPhiYdEta{element}(xi,eta).^2)./(g{element}(xi,eta).^2);
                g12{element} = @(xi,eta) (-dPhiXdEta{element}(xi,eta).*dPhiXdXi{element}(xi,eta) -dPhiYdEta{element}(xi,eta).*dPhiYdXi{element}(xi,eta))./(g{element}(xi,eta).^2);
                g22{element} = @(xi,eta) (dPhiXdXi{element}(xi,eta).^2 + dPhiYdXi{element}(xi,eta).^2)./(g{element}(xi,eta).^2);

                % Metric Tensor
                % g^{22}
%                 g22{element} = @(xi,eta) ((-mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 ...
%                                 + (mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2)./ ... 
%                                 (((-mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 ...
%                                 + (mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2).* ...
%                                 ((mapEta_Coeff1(elementEta)*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 + ...
%                                 (mapEta_Coeff1(elementEta)*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2));

                % g^{12} =  g^{21}
%                 g12{element} = @(xi,eta) zeros(size(xi));

                % g^{11}
%                 g11{element} = @(xi,eta) ((mapEta_Coeff1(elementEta)*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 + ...
%                                 (mapEta_Coeff1(elementEta)*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2)./ ... 
%                                 (((-mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 ...
%                                 + (mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2).* ...
%                                 ((mapEta_Coeff1(elementEta)*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 + ...
%                                 (mapEta_Coeff1(elementEta)*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2));

%                 g{element} = @(xi,eta) sqrt(((-mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 ...
%                                 + (mapXi_Coeff1(elementXi)*(mapEta_Coeff1(elementEta)*eta + mapEta_Coeff2(elementEta)).*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2).* ...
%                                 ((mapEta_Coeff1(elementEta)*cos(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2 + ...
%                                 (mapEta_Coeff1(elementEta)*sin(mapXi_Coeff1(elementXi)*xi + mapXi_Coeff2(elementXi))).^2));
            end
        end
        
    end

end