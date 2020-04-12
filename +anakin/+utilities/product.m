%{
product: a function to multiply two arrays and perform arbitrary
contractions. The resulting array indices are the free indices not
contracted, in the same order as given. 

T3 = anakin.utilities.product(T1,T2,<list1>,<list2>,...)

where: 
- <> denotes optional arguments
- T1,T2: the arrays to be multiplied
- list1,list2,...: vectors with the lists of subindices in T1*T2 (whole
  product) for each contraction.
  
AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
function T4 = product(T1,T2,varargin)

%% Compute sizes, ndims, numels
s1 = size(T1); % size, padded with ones to the right
s1_ = s1; % size, remove padding
for i = length(s1):-1:1
    if s1_(i) == 1
        s1_ = s1_(1:i-1);
    else
        break
    end
end
n1 = length(s1_); % actual number of dimensions 
e1 = numel(T1);
 
s2 = size(T2);
s2_ = s2;
for i = length(s2):-1:1
    if s2_(i) == 1
        s2_ = s2_(1:i-1);
    else
        break
    end
end
n2 = length(s2_); % actual number of dimensions 
e2 = numel(T2);

%% Whole product size, ndim
s3_ = [s1_,s2_]; % unpadded
s3 = [s3_,1,1]; % padded with 1s
n3 = n1+n2; 

%% Contraction mask
smask3 = true(1,n3); % mask of surviving dimensions
cmask3 = cell(1,nargin-2); % allocate for each contraction
for ic = 1:(nargin-2)
    cmask3{ic} = false(1,n3); 
    if all(varargin{ic}>0) && all(varargin{ic}<=n3) % skip bad contractions
        cmask3{ic}(varargin{ic}) = 1; 
    end
    smask3 = smask3 & ~cmask3{ic};
end
smask1 = smask3(1:n1);
smask2 = smask3(n1+1:end);

%% Contrated product size
s4_ = s3_(smask3); % unpadded
s4 = [s4_,1,1]; % padded with 1s 

%% Compute product 
T4 = zeros(s4); % allocate contracted product
if isa(T1,'sym') || isa(T2,'sym')
    T4 = sym(T4);
end
P3 = cell(1,n3); % subindices in the whole product
for i = 1:(e1*e2) % loops over the whole product
    [P3{:}] = ind2sub(s3,i); % get subindices from linear index
    skip = 0;
    for ic = 1:(nargin-2)
        cindices = [P3{cmask3{ic}}];
        if ~isempty(cindices) && ~all(cindices(end) == cindices) 
            skip = 1;
            break;
        end
    end 
    if skip
        continue;
    end
    P1 = P3(1:n1); % split P for T1 and T2
    P2 = P3(n1+1:end); 
    L1 = sub2ind(s1,P1{:},1,1); % get linear index from subindices
    L2 = sub2ind(s2,P2{:},1,1);    
    L4 = sub2ind(s4,P1{smask1},P2{smask2},1,1);
    T4(L4) = T4(L4) + T1(L1)*T2(L2);  
end










