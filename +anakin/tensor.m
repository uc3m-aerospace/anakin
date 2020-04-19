%{
DESCRIPTION:
tensor: class for creating and handling:
- scalars (0th-order tensors), 
- vectors (1st-order tensors), 
- 2nd-order square tensors, and 
- general nth-order square tensors. 
Many Matlab array functions have been overloaded for convenience. Some of
them only operate on tensors of a particular order, however.

SYNTAX:
T0 = anakin.tensor();  % returns default object 
T  = anakin.tensor(<T|A|c>,<B1>);
where: 
- <> denotes optional arguments
- | denotes alternative arguments
- T0 is the default tensor (scalar 0)
- T is a tensor  
- A is a point
- c is an array with the tensor components 
- B1 is a basis. If given, all previous input as relative to that basis
 
METHODS:
* spacedim: returns dimensionality of space
* components: returns the components of the tensor in a chosen basis
* x: returns a single component of the tensor in a chosen basis
* product: general product of tensors (with or without contraction)
* isscalar: true if ndims == 0
* isvector: true if ndims == 0
* is2tensor: true if ndims == 2
* dt: returns the time derivative of a tensor in a chosen basis
  (symbolic variables must be used)  
* subs: takes values of the symbolic unknowns and returns a tensor with
  purely numeric components (symbolic variables must be used)   

(for vectors only):
* magnitude: returns norm
* angle: returns angle between vectors
* unitvector: returns unit vector
* isunitvector: true if vector is unit vector
* isperpendicular: true if vectors are perpendicular
* isparallel: true if vectors are parallel
* plot: plots vector. Use FileExchange arrow3 code if available

(for 2nd order tensors only):
* ishermitian: true if tensor is hermitian
* isantihermitian: true if tensor is antihermitian
* isunitary: true if tensor is unitary
* isnormal: true if tensor is normal

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef tensor
    properties (Hidden = true, Access = protected) 
        c = 0; % components in the canonical vector basis
    end 
    methods % creation
        function T = tensor(c,B) % constructor
            % components
            if ~exist('c','var')
                c = 0; % default tensor
            elseif isa(c,'anakin.tensor') % tensor input
                c = c.c; 
            elseif isa(c,'anakin.point') % point input
                c = c.pos.c; 
            end 
            T.c = c;
            % change of basis
            if exist('B','var')   
                Binv = inv(B.matrix);
                for i = 1:T.ndims
                    T.c = anakin.utilities.product(T.c,Binv,[1,T.ndims+1]);
                end 
            end
        end
        function T = set.c(T,value) % on setting c
            % Simplify if needed
            if isa(value,'sym') % symbolic input
                value = formula(simplify(value)); % simplify and force sym rather than symfun to allow indexing into c
            end
            % Size of input
            sc = size(value);
            sc_ = sc(sc>1);       
            if ~isempty(sc_) && ~all(sc_(1) == sc_)
                error('The tensor class expects a square array as input');
            end
            % Assign
            T.c = reshape(value,[sc_,1,1]); % remove any singleton dimensions that may exist 
            try
                T.c = double(T.c); % turn sym without variables into double
            catch
                % pass
            end
        end
    end 
    methods (Hidden = true) % overloads
        function s = size(T) % size of non-singleton dimensions
            s = size(T.c); % size, padded with ones to the right
            for i = length(s):-1:1 % remove padding
                if s(i) == 1
                    s = s(1:i-1);
                else
                    break
                end
            end
        end
        function n = numel(T) % number of elements
            n = numel(T.c); 
        end
        function n = ndims(T) % number of non-singleton dimensions 
            n = length(size(T));
        end        
        function value = eq(T1,T2) % overload ==
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            if isa(T1.c,'sym') || isa(T2.c,'sym') % symbolic inputs
                value = isAlways(T1.c==T2.c,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                value = (abs(T1.c - T2.c) < 10*eps(T1.c)+10*eps(T2.c)); 
            end
            value = all(value(:));
        end
        function value = ne(T1,T2) % overload ~=
            value = ~eq(T1,T2);
        end
        function T3 = plus(T1,T2) % overloaded +
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = anakin.tensor(T1.c + T2.c);
        end
        function T3 = minus(T1,T2) % overloaded -
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = anakin.tensor(T1.c - T2.c);
        end 
        function T = uplus(T) % overloaded +(unitary)
            % pass
        end
        function T = uminus(T) % overloaded -(unitary)
            T.c = -T.c; 
        end 
        function T3 = mtimes(T1,T2) % overloaded *
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = product(T1,T2,T1.ndims+[0,1]); % contracts the innermost indices
        end 
        function T3 = mrdivide(T1,T2) % overloaded /
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = product(T1,inv(T2),T1.ndims+[0,1]);
        end
        function T3 = mldivide(T1,T2) % overloaded \
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = product(inv(T1),T2,T1.ndims+[0,1]);
        end
        function T = mpower(T,exponent) % overloaded ^
            T.c = T.c^exponent;
        end 
        function T = inv(T) % overloaded inv
            T.c = inv(T.c);
        end 
        function T = sin(T) % overloaded sin
            T.c = sin(T.c);
        end 
        function T = cos(T) % overloaded cos
            T.c = cos(T.c);
        end 
        function T = tan(T) % overloaded tan
            T.c = tan(T.c);
        end 
        function T = asin(T) % overloaded asin
            T.c = asin(T.c);
        end 
        function T = acos(T) % overloaded acos
            T.c = acos(T.c);
        end 
        function T = atan(T) % overloaded atan
            T.c = atan(T.c);
        end 
        function T3 = atan2(T1,T2) % overloaded atan2
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3.c = anakin.tensor(atan2(T1.c,T2.c));
        end 
        function T = exp(T) % overloaded exp
            T.c = exp(T.c);
        end 
        function T = log(T) % overloaded log
            T.c = log(T.c);
        end 
        function T = log10(T) % overloaded log10
            T.c = log10(T.c);
        end 
        function T = sign(T) % overloaded sign
            T.c = sign(T.c);
        end
        function T = abs(T) % overloaded abs
            T.c = abs(T.c);
        end 
        function T = sqrt(T) % overloaded sqrt
            T.c = sqrt(T.c);
        end
        function T = sum(T,varargin) % overloaded sum
            T.c = sum(T.c,varargin{:});
        end 
        function T = cumsum(T,varargin) % overloaded cumsum
            T.c = cumsum(T.c,varargin{:});
        end 
        function T = prod(T,varargin) % overloaded prod
            T.c = prod(T.c,varargin{:});
        end 
        function T = cumprod(T,varargin) % overloaded cumprod
            T.c = cumprod(T.c,varargin{:});
        end  
        function [lambdas,vectors] = eigs(T) % eigenvalues and eigenvectors
           [V,D] = eigs(T.c); 
           lambdas = cell(1,T.spacedim);
           vectors = cell(1,T.spacedim);
           for i = 1:T.spacedim
               lambdas{i} = anakin.tensor(D(i,i)); 
               vectors{i} = anakin.tensor(V(:,i));  
           end
        end
        function T = trace(T) % overloaded trace
            T.c = trace(T.c);
        end         
        function T = norm(T) % 2-norm
            T.c = norm(T.c);
        end
        function T3 = dot(T1,T2) % dot product
            T3 = T1*T2;
        end      
        function T3 = cross(T1,T2) % cross product
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            if T1.ndims ~= 1 || T2.ndims ~= 1 || T1.spacedim ~= 3 || T2.spacedim ~= 3
                error('This functionality is only available for vectors in 3D space');
            end
            T3 = anakin.tensor(cross(T1.c,T2.c)); 
        end
        function disp(T) % display
            switch T.ndims
                case 0
                    str = 'Scalar with value:';
                case 1
                    str = 'Vector with components:';
                case 2
                    str = 'Second-order tensor with components:';
                case 3
                    str = 'Third-order tensor with components:';
                otherwise
                    str = [num2str(T.ndims),'th-order tensor with components:'];
            end                    
            disp(str);
            disp(T.c);
        end
    end 
    methods % general functionality
        function n = spacedim(T) % dimensions of space (for first tensor dimension)
            s = [size(T),3]; % if size is empty (scalar case), default to 3.
            n = s(1); 
        end
        function c = components(T,B) % return components of T in basis B
            c = T.c; % if no basis is given, use the canonical vector basis
            if exist('B','var') 
                Bm = anakin.basis(B).matrix;
                for idim = 1:T.ndims
                    c = anakin.utilities.product(c,Bm,[1,T.ndims+1]);
                end
                if isa(c,'sym')
                    c = formula(simplify(c));
                end
            end
        end
        function x = x(T,varargin) % returns a single component x in basis B
            if isa(varargin{end},'anakin.basis') % last input may be B
                B = varargin{end};
                varargin = varargin(1:end-1);
            else
                B = anakin.basis(eye(T.spacedim)); % canonical vector basis
            end
            components = T.components(B); 
            x = components(varargin{:});                
        end
        function T3 = product(T1,T2,varargin) % Product of two ntensors, optionally with arbitrary contraction 
            if ~isa(T1,'anakin.tensor'); T1 = anakin.tensor(T1); end % ensure tensor
            if ~isa(T2,'anakin.tensor'); T2 = anakin.tensor(T2); end
            T3 = anakin.tensor(anakin.utilities.product(T1.c,T2.c,varargin{:})); 
        end
        function isscalar = isscalar(T)
            isscalar = (T.ndims == 0);
        end
        function isvector = isvector(T)
            isvector = (T.ndims == 1);
        end
        function is2tensor = is2tensor(T)
            is2tensor = (T.ndims == 2);
        end
        function dT = dt(T,B) % time derivative with respect to basis B. Requires symbolic tensor
            if exist('B','var')  
                dT = anakin.tensor(diff(sym(T.components(B)),1),B);
            else
                dT = anakin.tensor(diff(sym(T.components),1));
            end            
        end
        function T = subs(T,variables,values) % particularize symbolic tensor
            T.c = subs(T.c,variables,values);
            try
                T.c = double(T.c);
            catch
                % pass
            end
        end        
    end    
    methods % vector functionality
        function magnitude = magnitude(T) % alias for norm 2
            magnitude = norm(T); 
        end
        function angle = angle(T1,T2,T3) % angle between two vectors. A third one can be given to resolve sign
            if T1.ndims ~= 1 || T2.ndims ~= 1
                error('This functionality is only available for vectors');
            end
            angle = anakin.tensor(acos(dot(T1.c,T2.c)/(norm(T1.c)*norm(T2.c))));
            if exist('T3','var')
                if T1.spacedim ~= 3 || T2.spacedim ~= 3 || T3.spacedim ~= 3
                    error('This functionality is only available in 3D space');
                end
                angle.c = angle.c * sign(dot(cross(T1.c,T2.c),T3.c));
            end 
        end
        function unitvector = unitvector(T) % T divided by norm(T)
            unitvector = T/norm(T);
        end 
        function isunitvector = isunitvector(T) % T has unit norm
            if isa(T.c,'sym') % symbolic inputs
                isunitvector = isAlways(norm(T.c)==1,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isunitvector = (abs(norm(T.c)-1)<eps(1)); 
            end 
        end
        function isperpendicular = isperpendicular(T1,T2) % vectors are perpendicular
            if T1.ndims ~= 1 || T2.ndims ~= 1
                error('This functionality is only available for vectors');
            end
            if isa(T1.c,'sym') || isa(T2.c,'sym') % symbolic inputs
                isperpendicular = isAlways(dot(T1.c,T2.c)==0,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isperpendicular = (abs(dot(T1.c,T2.c))<eps(max(abs(T1.c(:))))+eps(max(abs(T2.c(:))))); 
            end 
        end
        function isparallel = isparallel(T1,T2) % vectors are parellel
            if T1.ndims ~= 1 || T2.ndims ~= 1
                error('This functionality is only available for vectors');
            end 
            if isa(T1.c,'sym') || isa(T2.c,'sym') % symbolic inputs
                isparallel = isAlways(cross(T1.c,T2.c) == [0;0;0],'Unknown','false'); % In case of doubt, false
            else
                isparallel = (cross(T1.c,T2.c) < eps(max(abs(T1.c(:))))+eps(max(abs(T2.c(:)))));
            end
            isparallel = all(isparallel);
        end
        function h = plot(T,varargin) % plot. First argument in varargin must be the O vector, if any
            if T.ndims ~= 1
                error('This functionality is only available for vectors');
            end 
            cc = [T.c; zeros(3-T.spacedim,1)];
            if mod(nargin,2) == 1 % no origin vector is given
                O = anakin.tensor([0;0;0]).c; % null vector
            else
                if isa(varargin{1},'anakin.point')
                    O = varargin{1}.pos.c;
                else
                    O = anakin.tensor(varargin{1}).c;
                end 
                varargin = varargin(2:end); % rest of varargin
            end 
            if exist('arrow3', 'file') % If Mathworks fileexchange arrow3 is available, use it
                hold on
                set(gca,'DataAspectRatio',[1,1,1]);
                h = arrow3(O',O'+cc',varargin{:});
                hold off
            else 
                h = line('XData',[O(1),O(1)+cc(1)],'YData',[O(2),O(2)+cc(2)],'ZData',[O(3),O(3)+cc(3)],'color','k'); 
                if ~isempty(varargin)
                    set(h,varargin{:}); % set options stored in varargin
                end
            end 
        end
    end
    methods % 2nd-order tensor functionality
        function ishermitian = ishermitian(T) % tensor is Hermitian (A' == A)
            if T.ndims ~= 2
                error('This functionality is only available for 2nd-order tensors');
            end
            if isa(T.c,'sym') % symbolic inputs
                ishermitian = isAlways(T.c' == T.c,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                ishermitian = (abs(T.c' - T.c) < eps(1)); 
            end 
            ishermitian = all(ishermitian(:));
        end  
        function isantihermitian = isantihermitian(T) % tensor is anti-Hermitian (A' == -A)
            if T.ndims ~= 2
                error('This functionality is only available for 2nd-order tensors');
            end
            if isa(T.c,'sym') % symbolic inputs
                isantihermitian = isAlways(T.c' == -T.c,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isantihermitian = (abs(T.c' + T.c) < eps(1)); 
            end 
            isantihermitian = all(isantihermitian(:));
        end  
        function isunitary = isunitary(T) % tensor is Unitary (A' * A = eye)
            if T.ndims ~= 2
                error('This functionality is only available for 2nd-order tensors');
            end
            if isa(T.c,'sym') % symbolic inputs
                isunitary = isAlways(T.c' * T.c == eye(T.spacedim,T.spacedim),'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isunitary = (abs(T.c' * T.c - eye(T.spacedim,T.spacedim)) < eps(1)); 
            end 
            isunitary = all(isunitary(:));
        end  
        function isnormal= isnormal(T) % tensor is Normal (A' * A = A * A')
            if T.ndims ~= 2
                error('This functionality is only available for 2nd-order tensors');
            end
            if isa(T.c,'sym') % symbolic inputs
                isnormal = isAlways(T.c' * T.c == T.c * T.c','Unknown','false'); % In case of doubt, false
            else % numeric input            
                isnormal = (abs(T.c' * T.c - T.c * T.c') < eps(1)); 
            end 
            isnormal = all(isnormal(:));
        end  
    end 
end








