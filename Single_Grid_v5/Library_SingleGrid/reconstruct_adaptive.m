function varargout = reconstruct_adaptive(form,varargin)

global Nadaptive numElements nn Nmax

if isempty(numElements); numElements = 1; end

switch form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 0  % Reconstruction zero-forms

        % uu = reconstruct_zeroforms(U,h)
        U = varargin{1};
        h = varargin{2};

        hh = cell(1,Nmax);
        for i=unique(Nadaptive)
            hh{i} = kron(h{i},h{i})';
        end

        uu = zeros(nn^2,numElements);
        for i=1:numElements
            uu(:,i) = hh{Nadaptive(i)}*U(1:(Nadaptive(i)+1)^2,i);
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

        he = cell(1,Nmax);
        for i=unique(Nadaptive)
            he{i} = kron(e{i},h{i})';
        end
        
        eh = cell(1,Nmax);
        for N=unique(Nadaptive)
        
            eh{N} = zeros(nn^2,N*(N+1));
            for i=1:N
                for j=1:N+1
                    ij = j+(i-1)*(N+1);
                    eh{N}(:,ij) = kron(h{N}(j,:),e{N}(i,:))';
                end
            end

        end

        qx = zeros(nn^2,numElements);
        qy = zeros(nn^2,numElements);
        for i=1:numElements

            N = Nadaptive(i);

            qxi  = he{N}*Qxi(1:N*(N+1),i);
            qeta = eh{N}*Qeta(1:N*(N+1),i);

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
        
        ee = cell(1,Nmax);
        for i=unique(Nadaptive)
            ee{i} = kron(e{i},e{i})';
        end
        
        uu = zeros(nn^2,numElements);
        for i=1:numElements
            uu(:,i) = (ee{Nadaptive(i)}*U(1:Nadaptive(i)^2,i))./Mesh.J(:,i);
        end

        varargout{1} = uu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end