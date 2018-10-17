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
function tests = test_point
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    S1 = frame([0,1,0;-1,0,0;0,0,1]); 
    
    A = point; 
    A = point(A); 
    A = point([1;2;0]); % coordinates in canonical reference frame    
    A = point([1,2,0]); % coordinates in canonical tensor basis, row array form    
    A = point(tensor([1;2;0])); % tensor 
    A = point([1;2;0],S1); % coordinates in another basis, column array form
    A = point([1,2,0],S1); % coordinates in another basis, row array form    
    A = point(tensor([1;2;0]),S1); % tensor     
end

function test_coordinates(~) % call coordinates, x, y, z
    import anakin.*
    S1 = frame([2,2,2],[0,1,0;-1,0,0;0,0,1]); 
    
    A = point([1,2,3],S1); 
    assert(all(A.coordinates(S1) == [1;2;3]));
    assert(A.x(1,S1) == 1);
    assert(A.x(2,S1) == 2);
    assert(A.x(3,S1) == 3);
end

function test_posvelaccel(~) % Call position, velocity and acceleration wrt canonical reference frame
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);
        B0 = basis;
        S0 = frame; 
        o = tensor([1,t,cos(theta)]);
        B = B0.rotatex(theta);
        S1 = frame(point(o),B);

        a = tensor([t,sin(theta),xi],B);
        A = point(o+a);
        
        assert(A.pos(S1) == tensor([t;sin(theta);xi],S1.basis));
        assert(A.vel(S1) == tensor([1;cos(theta)*diff(theta,t);diff(xi)],S1.basis));
        assert(A.accel(S1) == tensor([0; cos(theta(t))*diff(theta(t), 2) - sin(theta(t))*diff(theta(t), t)^2; diff(xi(t), 2)],S1.basis));

        assert(A.pos(S0) == o + a); 
        assert(A.vel(S0) == o.dt(S0.basis) + a.dt(S0.basis));
        assert(A.accel(S0) == o.dt(S0.basis).dt(S0.basis) + a.dt(S0.basis).dt(S0.basis));
    end
end 

function test_subs(~) % Particularize a symbolic tensor
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) phi(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'),in(xi(t), 'real')]);    
        A = point([cos(theta),sin(theta)+xi^2*sin(phi),xi*cos(phi)]); 

        B = A.subs({t,theta,phi,xi},{1,t^2-4,2*t+3,-2}); % call with cell arrays
        B = A.subs([t,theta,phi,xi],[1,t^2-4,2*t+3,-2]); % call with arrays
    end
end
  
function test_plot(~) % tensor plotting 
    import anakin.*
    A = point([1,5,3]); 
    
    f = figure;
    A.plot;
    A.plot('color','b');
    A.plot('color',[0.5,0.5,0.5]);
    close(f);
end







