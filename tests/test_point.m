%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. For convenience, you can set the run
  botton of Matlab's interface to runtests(<this function>)
* Your working directory must be the
  directory where this test file is contained.  
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
    
    P = point;
end

function test_creator2(~) % Call creator with numeric arguments
    import anakin.*
    S1 = frame([0,1,0;-1,0,0;0,0,1]); 
    
    P = point([1;2;0]); % coordinates in canonical reference frame
    P = point([1,2,0]); % coordinates in canonical vector basis, row array form    
    P = point([1;2;0],S1); % coordinates in another basis, column array form
    P = point([1,2,0],S1); % coordinates in another basis, row array form    
    P = point(1,2,0); % coordinates in canonical vector basis, independently given
    P = point(1,2,0,S1); % coordinates in another basis, independently given
end

function test_posvelaccel(~) % Call position, velocity and acceleration wrt canonical reference frame
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);
        B0 = basis;
        S0 = frame; 
        o = vector(1,t,cos(theta));
        B = B0.rotatex(theta);
        S1 = frame(point(o),B);

        p = vector(t,sin(theta),xi,B);
        P = point(o+p);

        assert(P.pos(S1) == vector([ t;sin(theta);xi],S1));
        assert(P.vel(S1) == vector([ 1;cos(theta)*diff(theta,t);diff(xi)],S1));
        assert(P.accel(S1) == vector([ 0; cos(theta(t))*diff(theta(t), 2) - sin(theta(t))*diff(theta(t), t)^2; diff(xi(t), 2)],S1));

        assert(P.pos(S0) == o + p); 
        assert(P.vel(S0) == o.dt(S0) + p.dt(S0));
        assert(P.accel(S0) == o.dt(S0).dt(S0) + p.dt(S0).dt(S0));
    end
end 

function test_subs(~) % Particularize a symbolic vector
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) phi(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'),in(xi(t), 'real')]);    
        P = point([cos(theta),sin(theta)+xi^2*sin(phi),xi*cos(phi)]); 

        Q = P.subs({t,theta,phi,xi},{1,t^2-4,2*t+3,-2}); % call with cell arrays
        Q = P.subs([t,theta,phi,xi],[1,t^2-4,2*t+3,-2]); % call with arrays
    end
end
  
function test_plot(~) % vector plotting 
    import anakin.*
    P = point([1,5,3]); 
    
    f = figure;
    P.plot;
    P.plot('color','b');
    P.plot('color',[0.5,0.5,0.5]);
    close(f);
end







