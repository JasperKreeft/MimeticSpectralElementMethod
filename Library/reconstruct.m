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
                eh(:,ij) = kron(h(j,:),e(i,:));
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
        if N==1 && size(U,1)>size(U,2); U = U'; end
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


    case 11 % Reconstruction covector-valued one-forms

        Txixi   = varargin{1};
        Tetaxi  = varargin{2};
        Txieta  = varargin{3};
        Tetaeta = varargin{4};
        hgl    = varargin{5};
        egl    = varargin{6};
        heg    = varargin{7};
        eeg    = varargin{8};        
        Mesh   = varargin{9};

        nn = size(hgl,2);

        he = kron(egl,heg)';
        eh = zeros(nn^2,(N+1)^2);
        for i=1:N+1
            for j=1:N+1
                ij = i+(j-1)*(N+1);
                eh(:,ij) = kron(hgl(j,:),eeg(i,:));
            end
        end

        txx = zeros(nn^2,numElements);
        tyx = zeros(nn^2,numElements);
        for i=1:numElements

            txixi  = he*Txixi(:,i);
            tetaxi = eh*Tetaxi(:,i);

            txx(:,i) = ( txixi.*Mesh.dXdXi(:,i) + tetaxi.*Mesh.dXdEta(:,i) )./Mesh.J(:,i);
            tyx(:,i) = ( txixi.*Mesh.dYdXi(:,i) + tetaxi.*Mesh.dYdEta(:,i) )./Mesh.J(:,i);


        end

        varargout{1} = txx;
        varargout{2} = tyx;

        he = repmat(kron(hgl',ones(1,N+1)),nn,1).*repmat(kron(eeg',ones(nn,1)),1,N+1);
        eh = zeros(nn^2,N*(N+1));
        for i=1:N
            for j=1:N+2
                ij = j+(i-1)*(N+2);
                eh(:,ij) = kron(heg(j,:),egl(i,:));
            end
        end

        txy = zeros(nn^2,numElements);
        tyy = zeros(nn^2,numElements);
        for i=1:numElements

            txieta  = he*Txieta(:,i);
            tetaeta = eh*Tetaeta(:,i);

            txy(:,i) = ( txieta.*Mesh.dXdXi(:,i) + tetaeta.*Mesh.dXdEta(:,i) )./Mesh.J(:,i);
            tyy(:,i) = ( txieta.*Mesh.dYdXi(:,i) + tetaeta.*Mesh.dYdEta(:,i) )./Mesh.J(:,i);


        end

        varargout{3} = txy;
        varargout{4} = tyy;        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 21 % Reconstruction covector-valued two-forms

        Mx = varargin{1};
        My = varargin{2};
        
        if N==1 && size(Mx,1)>size(Mx,2); Mx = Mx'; My = My'; end
        
        egl = varargin{3};
        eeg = varargin{4};
        Mesh = varargin{5};

        nn = size(egl,2);

        eegegl = kron(egl,eeg)'; 
        egleeg = zeros(nn^2,N*(N+1));
        for i=1:N
            for j=1:N+1
                ij = j+(i-1)*(N+1);
                egleeg(:,ij) = kron(eeg(j,:),egl(i,:));
            end
        end

        mmx = zeros(nn^2,numElements);
        mmy = zeros(nn^2,numElements);
        for i=1:numElements
            mmx(:,i) = (eegegl*Mx(:,i))./Mesh.J(:,i);
            mmy(:,i) = (egleeg*My(:,i))./Mesh.J(:,i);
        end

        varargout{1} = mmx;
        varargout{2} = mmy;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end