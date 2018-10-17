%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. For convenience, you can set the run
  botton of Matlab's interface to runtests(<this function>)
* Your working directory must be the
  directory where this test file is contained.  

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
function tests = test_product
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_whole(~) 
    product = @anakin.utilities.product;
    T1 = 2;
    T3 = [1;2;3];
    T13 = [1,2,3];
    T33 = [1,2,3;4,5,6;7,8,9]; 
     
    assert(product(T1,T1) == 4); % two scalars
    assert(all(product(T1,T3) == T1*T3)); % scalar-vector
    assert(all(product(T3,T1) == T1*T3)); % vector-scalar
    assert(all(all(product(T3,T3) == T3*T3'))); % vector-vector
    assert(all(product(T1,T13) == T1*T13)); % scalar-vector
    assert(all(product(T13,T1) == T1*T13)); % vector-scalar
    assert(all(all(product(T1,T33) == T1*T33))); % scalar-matrix
    assert(all(all(product(T33,T1) == T33*T1))); % matrix-scalar
    
    % vector-vector
    T = product(T13,T13);
    Texpected(1,:,1,:) = T13'*T13;
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
        
    % vector-vector
    T = product(T3,T13);
    Texpected(:,1,:) = T3*T13;
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
    
    % vector-vector
    T = product(T13,T3);
    Texpected(1,:,:) = T13'*T3';
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
    
    % matrix-vector
    T = product(T33,T3);
    Texpected(:,:,1) = T33*T3(1); Texpected(:,:,2) = T33*T3(2); Texpected(:,:,3) = T33*T3(3);
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
    
    % vector-matrix
    T = product(T3,T33);
    Texpected(1,:,:) = T33*T3(1); Texpected(2,:,:) = T33*T3(2); Texpected(3,:,:) = T33*T3(3);
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
    
    % matrix-matrix
    T = product(T33,T33);
    Texpected(:,:,1,1) = T33*T33(1,1); Texpected(:,:,1,2) = T33*T33(1,2); Texpected(:,:,1,3) = T33*T33(1,3); 
    Texpected(:,:,2,1) = T33*T33(2,1); Texpected(:,:,2,2) = T33*T33(2,2); Texpected(:,:,2,3) = T33*T33(2,3); 
    Texpected(:,:,3,1) = T33*T33(3,1); Texpected(:,:,3,2) = T33*T33(3,2); Texpected(:,:,3,3) = T33*T33(3,3); 
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
end

function test_contraction(~) 
    product = @anakin.utilities.product;
    T3 = [1;2;3];
    T13 = [1,2,3];
    T33 = [1,2,3;
           4,5,6; 
           7,8,9]; 
     
    assert(product(T3,T3,[1,2]) == 14); % two vectors 
    assert(product(T3,T13,[1,3]) == 14); % two vectors
    assert(product(T13,T3,[2,3]) == 14); % two vectors
    assert(product(T13,T13,[2,4]) == 14); % two vectors
    assert(all(all((product(T13,T13,[1,3]) == T13'*T13)))); % two vectors
    
    % vector-matrix
    T = product(T3,T33,[1,2]);
    Texpected(1,:) = [30,36,42];
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
    
    % vector-matrix
    T = product(T3,T33,[1,3]);
    Texpected(:) = [14;32;50];
    assert(all(T(:) == Texpected(:))); 
    clear Texpected
   
    % matrix-matrix
    T = product(T33,T33,[1,3]);
    Texpected = T33'*T33;
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
    
    % matrix-matrix
    T = product(T33,T33,[2,3]);
    Texpected = T33*T33;
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
    
    % matrix-matrix
    T = product(T33,T33,[1,3],[2,4]);
    Texpected = 285;
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
    
    % matrix-matrix
    T = product(T33,T33,[1,4],[2,3]);
    Texpected = 261;
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
    
    % matrix self contract
    T = product(T33,1,[1,2]);
    Texpected = 15;
    assert(all(T(:) == Texpected(:)));    
    clear Texpected
    
end





















