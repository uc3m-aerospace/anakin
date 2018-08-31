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
function tests = test_particle
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    S1 = frame([4,5,6],[0,1,0;-1,0,0;0,0,1]); 
    x = 1; y = 2; z = 3; c = [x;y;z]; cp = [6;4;9];
    a = vector(c); A = point(c); 
    mass = 3;
    
    P = particle();
    assert(all(P.pos.components == [0;0;0]));
    assert(all(P.mass == 1));
    
    P = particle(particle);
    assert(all(P.pos.components == [0;0;0]));
    assert(all(P.mass == 1));
    
    P = particle(A);
    assert(all(P.pos.components == c));
    assert(all(P.mass == 1));
    
    P = particle(a);
    assert(all(P.pos.components == c));
    assert(all(P.mass == 1));
    
    P = particle(c);
    assert(all(P.pos.components == c));
    assert(all(P.mass == 1));
    
    P = particle(c');
    assert(all(P.pos.components == c));
    assert(all(P.mass == 1));
    
    P = particle(x,y,z);
    assert(all(P.pos.components == c));
    assert(all(P.mass == 1));
    
    P = particle(particle,S1);
    assert(all(P.pos.components == [4;5;6]));
    assert(all(P.mass == 1));
    
    P = particle(A,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == 1));
    
    P = particle(a,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == 1));
    
    P = particle(c,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == 1));
    
    P = particle(c',S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == 1));
    
    P = particle(x,y,z,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == 1));
     
    P = particle(mass);
    assert(all(P.pos.components == [0;0;0]));
    assert(all(P.mass == mass));
    
    P = particle(mass,A);
    assert(all(P.pos.components == c));
    assert(all(P.mass == mass));
    
    P = particle(mass,a);
    assert(all(P.pos.components == c));
    assert(all(P.mass == mass));
    
    P = particle(mass,c);
    assert(all(P.pos.components == c));
    assert(all(P.mass == mass));
    
    P = particle(mass,c');
    assert(all(P.pos.components == c));
    assert(all(P.mass == mass));
    
    P = particle(mass,x,y,z);
    assert(all(P.pos.components == c));
    assert(all(P.mass == mass));
     
    P = particle(mass,A,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == mass));
    
    P = particle(mass,a,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == mass));
    
    P = particle(mass,c,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == mass));
    
    P = particle(mass,c',S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == mass));
    
    P = particle(mass,x,y,z,S1);
    assert(all(P.pos.components == cp));
    assert(all(P.mass == mass));
    
end

function test_momentum(~) % call momentum, linear momentum, energy
    import anakin.*
    if license('test','symbolic_toolbox')
        syms t;
        syms theta(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real'),...
                               in(diff(theta(t),t), 'real'), in(diff(xi(t),t), 'real')]);    
        c = formula([cos(theta);xi^2;xi]); 
        cp = formula([-sin(theta)*diff(theta,t);2*xi*diff(xi,t);diff(xi,t)]); 
        mass = 3;    
        P = particle(mass,c);

        assert(all(isAlways(P.p.components == 3*cp)));
        assert(all(isAlways(P.H.components == 3*cross(c,cp))));
        assert(all(isAlways(P.T == 3*dot(cp,cp)/2)));
    end  
end

function test_inertia(~) % Inertia
    import anakin.*
    S1 = frame([4,5,6],[0,1,0;-1,0,0;0,0,1]); 
    x = 1; y = 2; z = 3; c = [x;y;z]; cp = [6;4;9]; 
    mass = 3;

    P = particle(mass,c);
    i = P.inertia(S1);
 
end

function test_subs(~) % Particularize a symbolic vector
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real'),...
                               in(diff(theta(t),t), 'real'), in(diff(xi(t),t), 'real')]);    
        c = formula([cos(theta);xi^2;xi]);  
        mass = 3;    
        P = particle(mass,c);

        assert(P.subs({t,theta,xi},{1,3,-2}).pos == vector([cos(3);4;-2])); % call with cell arrays 
    end
end