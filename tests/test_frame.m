%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. For convenience, you can set the run
  botton of Matlab's interface to runtests(<this function>)
* Your working directory must be the
  directory where this test file is contained.  
%}
function tests = test_frame
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    
    S = frame;
end

function test_creator2(~) % Call creator with arguments
    import anakin.*
    syms t;
    syms theta(t) xi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);
    O = vector([xi(t),0,0]);
    i = vector([cos(theta),sin(theta),0]);
    j = vector([-sin(theta),cos(theta),0]);
    k = vector([0,0,1]);
    matrix = formula([cos(theta),-sin(theta),0; sin(theta),cos(theta),0; 0,0,1]);
    B = basis(i,j,k); 
    S1 = frame;
    
    S = frame(frame); % repetition    
    S = frame(O); % origin         
    S = frame(B); % basis
    S = frame([1,2,3]); % origin components 
    
    S = frame(O,B); % origin, basis. THIS IS THE RECOMMENDED WAY TO CREATE A FRAME
    S = frame([1,2,3],S1); % relative components of origin and frame 
    
    S = frame(1,2,3); % x y z
    S = frame(i,j,k); % i j k
    S = frame([1,2,3],B,S1); % relative origin components, basis, frame  
    
    S = frame(O,i,j,k); % origin, ijk
    S = frame([1,2,3],i,j,k); % origin components, ijk
    S = frame(1,2,3,B); % xyz, basis 
    S = frame(1,2,3,S1); % xyz, frame
    S = frame(i,j,k,S1); % ijk, frame
    
    S = frame(1,2,3,B,S1); % xyz, basis, frame 
    S = frame([1,2,3],i,j,k); % origin components, ijk, frame 
    
    S = frame(1,2,3,i,j,k,S1); % xyz, ijk, frame        
end

function test_overloads(~) % overloaded operators
    import anakin.*
    syms t;
    syms theta(t) phi(t) xi(t) eta(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'), in(xi(t), 'real'), in(eta(t), 'real')]);    
    O1 = vector([xi(t),0,0]);
    O2 = vector([eta(t),xi(t)*eta(t),0]) + O1; % defined wrt O1
    B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]); 
    B2 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1],B1);     
    S1 = frame(O1,B1);
    S2 = frame(O2,B2);

    assert(S1==S1);
    assert(S1~=S2);
end

function test_posvelaccel(~) % Call position, velocity and acceleration wrt canonical reference frame
    import anakin.*
    syms t;
    syms theta(t) xi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);    
    O = point([xi(t),0,0]); 
    B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);         
    S = frame(O,B); 
    
    assert(S.pos == vector([ xi;0;0]));
    assert(S.vel == vector([ diff(xi(t), 1);0;0]));
    assert(S.accel == vector([ diff(xi(t), 2);0;0]));
end

function test_omegaalpha(~) % Call angular velocity and acceleration wrt canonical reference frame
    import anakin.*
    syms t;
    syms theta(t) xi(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);    
    O = vector([xi(t),0,0]);
    B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);
    S = frame(O,B); 
    
    assert(S.omega == vector([ 0;0; diff(theta(t), 1)]));
    assert(S.alpha == vector([ 0;0; diff(theta(t), 2)]));
end

function test_omegaalpha2(~) % Call angular velocity, angular acceleration wrt a second reference frame
    import anakin.*
    syms t;
    syms theta(t) phi(t) xi(t) eta(t);
    assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'), in(xi(t), 'real'), in(eta(t), 'real')]);    
    O1 = vector([xi(t),0,0]);
    O2 = vector([eta(t),xi(t)*eta(t),0]) + O1; % defined wrt O1
    B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]); 
    B2 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1],B1); % defined wrt B1
    S1 = frame(O1,B1);
    S2 = frame(O2,B2);
     
    assert(S2.omega(S1) == vector([ 0;0; diff(theta(t), 1)],B1));
    assert(S2.alpha(S1) == vector([ 0;0; diff(theta(t), 2)],B1));
end

function test_origin(~) % origin point
    import anakin.*
    S = frame(1,2,3);
    
    origin = S.origin;
end

function test_subs(~) % subs
    import anakin.*
    syms t;
    syms phi(t) xi(t);
    assume([in(t, 'real'), in(phi(t), 'real'), in(xi(t), 'real')]);    
    O1 = vector([xi(t),0,0]); 
    B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);
    S1 = frame(O1,B1); 
     
    S1_ = S1.subs({t,phi,xi},{2,2*t,3*t});
end

function test_rotated_displaced(~) % create rotated and displaced frames
    import anakin.*
    syms t;
    syms mypsi(t) theta(t) phi(t) xi(t);
    assume([in(t, 'real'), in(mypsi(t), 'real'), in(theta(t), 'real'), in(phi(t), 'real'), in(xi(t), 'real')]);
    S0 = frame;
    
    Sx = S0.rotatex(pi/6); % with numeric value
    Sx = S0.rotatex(theta); % with symbolic expression
    Sy = S0.rotatey(theta);
    Sz = S0.rotatez(theta);
 
    S1 = S0.rotatez(mypsi).rotatex(theta).rotatez(phi);    
    
    Sd = S0.displace(vector(1,2,3));
end

 
function test_plot(~) % vector plotting  
    import anakin.*
    S = frame;
    
    f = figure;
    S.plot; 
    close(f);
end











