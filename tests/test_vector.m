%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. For convenience, you can set the run
  botton of Matlab's interface to runtests(<this function>)
* Your working directory must be the
  directory where this test file is contained.  
%}
function tests = test_vector
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    
    a = vector;
end

function test_creator2(~) % Call creator with numeric arguments
    import anakin.*
    B1 = basis([0,1,0;-1,0,0;0,0,1]); 
    
    a = vector([1;2;0]); % components in canonical vector basis, column array form
    a = vector([1,2,0]); % components in canonical vector basis, row array form    
    a = vector([1;2;0],B1); % components in another basis, column array form
    a = vector([1,2,0],B1); % components in another basis, row array form    
    a = vector(1,2,0); % components in canonical vector basis, independently given
    a = vector(1,2,0,B1); % components in another basis, independently given
end

function test_creator3(~) % Call creator with arguments of type sym
    import anakin.*
    syms t;
    syms theta(t) phi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
    B1 = basis([1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)]);
    
    a = vector([cos(theta);sin(theta);0]); % components in canonical vector basis, column array form
    a = vector([cos(theta),sin(theta),0]); % components in canonical vector basis, row array form
    a = vector([cos(theta);sin(theta);0],B1); % components in another basis, column array form
    a = vector([cos(theta),sin(theta),0],B1); % components in another basis, row array form
    a = vector(cos(theta),sin(theta),0); % components in canonical vector basis, independently given
    a = vector(cos(theta),sin(theta),0,B1); % components in another basis, independently given
end

function test_components(~) % Call components, x,y,z
    import anakin.*
    syms t;
    syms theta(t) phi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
    B1 = basis([1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)]);    
    a = vector([cos(theta),sin(theta),0],B1);  
    
    assert(isAlways(all(a.components(B1) == [cos(theta);sin(theta);0])));
    assert(isAlways(a.x(B1) == cos(theta)));
    assert(isAlways(a.y(B1) ==  sin(theta)));
    assert(isAlways(a.z(B1) == 0));
end

function test_overloads(~) % components, x,y,z
    import anakin.*
    syms t;
    syms theta(t) phi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);    
    a = vector([cos(theta),sin(theta),0]); 
    b = vector([0,cos(phi),sin(phi)]);     
    
    a+b;a-b;+a;-a;2*a;a*2; 2.1 .*a;a.*2;a/2;2\a;a./2;21.2 .\a; % many overloaded operators
    dot(a,b);
    norm(a);
    cross(a,b);
    dir(a);
    magnitude(a);  
    assert(a==a);
    assert(a~=b);
end

function test_isa(~) % isunitary, isparallel, isperpendicular
    import anakin.*
    syms t;
    syms theta(t) phi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);    
    a = vector([cos(theta),sin(theta),0]); 
    b = vector([0,0,sin(phi)]);         
    c = 2*a; 
    
    assert(isunitary(a)); 
    assert(~isunitary(c)); 
    assert(isperpendicular(a,b));
    assert(~isperpendicular(a,c));
    assert(~isparallel(a,b));
    assert(isparallel(a,c)); 
end

function test_subs(~) % Particularize a symbolic vector
    import anakin.*
    syms t;
    syms theta(t) phi(t) xi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'),in(xi(t), 'real')]);    
    a = vector([cos(theta),sin(theta)+xi^2*sin(phi),xi*cos(phi)]); 
    
    c = a.subs({t,theta,phi,xi},{1,t^2-4,2*t+3,-2}); % call with cell arrays
    c = a.subs([t,theta,phi,xi],[1,t^2-4,2*t+3,-2]); % call with arrays
end

function test_timediff(~) % Time derivative wrt a given basis calling dt
    import anakin.*
    syms t;
    syms theta(t) phi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
    B0 = basis;
    B1 = basis([1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)]);    
    a = vector([cos(theta),sin(theta),0],B1); 
    
    assert(a.dt(B1) == vector(diff(theta,1)*[-sin(theta),cos(theta),0],B1)); % time derivative of a wrt B1
    assert(a.dt(B0) == a.dt(B1) + cross(omega(B1,B0),a)); % time derivative of a wrt B0
end

function test_plot(~) % vector plotting 
    import anakin.*
    a = vector([1,5,3]); 
    
    f = figure;
    a.plot;
    a.plot('color','b');
    a.plot(vector([1,1,1]),'color',[0.5,0.5,0.5]);
    close(f);
end







