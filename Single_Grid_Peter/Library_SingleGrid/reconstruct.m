function varargout = reconstruct(form,varargin)

global N numElements

if isempty(numElements); numElements = 1; end

switch form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 0  % Reconstruction zero-forms

        % uu = reconstruct_zeroforms(U,h)
        U = varargin{1};
        h = varargin{2};

        nn = size(h,2);

        hh = kron(h,h)';

        uu = zeros(nn^2,numElements);
        for i=1:numElements
            uu(:,i) = hh*U(:,i);
        end
        
        varargout{1} = uu;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 1  % Reconstruction one-forms
        
        % function [qx,qy,qMag] = reconstruct_oneforms(Qxi,Qeta,h,e,Mesh)
        Qxi  = varargin{1};
        Qeta = varargin{2};
        h    = varargin{3};
        e    = varargin{4};
        Mesh = varargin{5};

        nn = size(h,2);

        he = kron(e,h)';
        eh = zeros(nn^2,N*(N+1));
        for i=1:N
            for j=1:N+1
                ij = j+(i-1)*(N+1);
                eh(:,ij) = kron(h(j,:),e(i,:))';
            end
        end

        qx = zeros(nn^2,numElements);
        qy = zeros(nn^2,numElements);
        for i=1:numElements

            qxi  = he*Qxi(:,i);
            qeta = eh*Qeta(:,i);

            qx(:,i) = ( qxi.*Mesh.dXdXi(:,i) + qeta.*Mesh.dXdEta(:,i) )./Mesh.J(:,i);
            qy(:,i) = ( qxi.*Mesh.dYdXi(:,i) + qeta.*Mesh.dYdEta(:,i) )./Mesh.J(:,i);


        end

        varargout{1} = qx;
        varargout{2} = qy;

        if nargout==3
            % qMag is magnitude of q
            qMag = sqrt(qx.^2+qy.^2);
            varargout{3} = qMag;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2  % Reconstruction two-forms
        
        % function uu = reconstruct_twoforms(U,e,J)
        U = varargin{1};
        e = varargin{2};
        Mesh = varargin{3};
        
        nn = size(e,2);

        ee = kron(e,e)';

        uu = zeros(nn^2,numElements);
        for i=1:numElements
            uu(:,i) = (ee*U(:,i))./Mesh.J(:,i);
        end

        varargout{1} = uu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end