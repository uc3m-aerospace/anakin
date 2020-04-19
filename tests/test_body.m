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
function tests = test_body
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*

    b = body();
    b = body(13,2*eye(3),[1;2;3],eye(3));
     
    b = body();
    b = body(b); % repetition
    for i = 1:2
        if i == 1
            t = 3.1;
            theta = pi/6;
            xi = 1.2;
        elseif license('test','symbolic_toolbox') 
            syms t;
            syms theta(t) xi(t);
            assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real')]);
        else 
            break;
        end 
  
        i = tensor([cos(theta),sin(theta),0]);
        j = tensor([-sin(theta),cos(theta),0]);
        k = tensor([0,0,1]); 
        B = basis(i,j,k); 
        S = frame;
        S1 = frame(B,[4,4,4]);
        M = tensor(22);
        P = point([xi,0,0]);
        
        b = body(P,M,S,B);
        b = body(P,M,S1);
        b = body(P,M,B,S);
        b = body(S1,S);
        b = body(S,S1);
        b = body(b,S);
        
        b = body(eye(3)*8);
        b = body(tensor(eye(3))*8);
        
    end
end
    
function test_momentum(~) % call momentum, linear momentum, energy
    import anakin.*
    if license('test','symbolic_toolbox')
        syms t;
        syms theta(t) xi(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real'), in(phi(t), 'real'),...
                in(diff(theta(t),t), 'real'), in(diff(xi(t),t), 'real'), in(diff(phi(t),t), 'real')]);    
        c = formula([cos(theta);xi^2;xi]); 
        cp = formula([-sin(theta)*diff(theta,t);2*xi*diff(xi,t);diff(xi,t)]); 
        a = tensor(c);
        mass = 3;    
        
        i = tensor([cos(phi),sin(phi),0]);
        j = tensor([-sin(phi),cos(phi),0]);
        k = tensor([0,0,1]); 
        B = basis(i,j,k); 
        
        I0 = rand(3); I0 = I0 + I0' + eye(3) * 8; I0 = tensor(I0,B);
        
        b = body(mass,a,B,I0);

        assert(all(isAlways(b.p.components == mass*cp)));
        assert(all(isAlways(b.H.components == mass*cross(c,cp) + I0.components * b.omega.components)));
        assert(all(isAlways(b.T == mass*dot(cp,cp)/2 + b.omega.components' * I0.components * b.omega.components / 2)));
    end  
end 

function test_subs(~) % Particularize a symbolic vector
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) xi(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real'), in(phi(t), 'real'),...
                in(diff(theta(t),t), 'real'), in(diff(xi(t),t), 'real'), in(diff(phi(t),t), 'real')]);    
        c = formula([cos(theta);xi^2;xi]); 
        cp = formula([-sin(theta)*diff(theta,t);2*xi*diff(xi,t);diff(xi,t)]); 
        a = tensor(c);
        mass = 3;    
        
        i = tensor([cos(phi),sin(phi),0]);
        j = tensor([-sin(phi),cos(phi),0]);
        k = tensor([0,0,1]); 
        B = basis(i,j,k); 
        
        I0 = rand(3); I0 = I0 + I0' + eye(3) * 8; I0 = tensor(I0,B);
        
        b = body(mass,a,B,I0);

        assert(b.subs({t,theta,xi},{1,3,-2}).pos == tensor([cos(3);4;-2])); % call with cell arrays 
    end
end

