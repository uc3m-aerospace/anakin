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
    S = frame(S); % repetition
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
 

        O = tensor([xi,0,0]);
        i = tensor([cos(theta),sin(theta),0]);
        j = tensor([-sin(theta),cos(theta),0]);
        k = tensor([0,0,1]); 
        B = basis(i,j,k); 
        S1 = frame;

        S = frame(frame); % repetition    
        S = frame(O); % origin         
        S = frame(B); % basis
        S = frame([1,2,3]); % origin components 

        S = frame(O,B); % origin, basis. THIS IS THE RECOMMENDED WAY TO CREATE A FRAME
        S = frame([1,2,3],S1); % relative components of origin and frame 
 
        S = frame([1,2,3],B,S1); % relative origin components, basis, frame   
    end
end

function test_overloads(~) % overloaded operators
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) phi(t) xi(t) eta(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'), in(xi(t), 'real'), in(eta(t), 'real')]);    
        O1 = tensor([xi(t),0,0]);
        O2 = tensor([eta(t),xi(t)*eta(t),0]) + O1; % defined wrt O1
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]); 
        B2 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1],B1);     
        S1 = frame(O1,B1);
        S2 = frame(O2,B2);

        assert(S1==S1);
        assert(S1~=S2);
    end
end
  
function test_subs(~) % subs
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms phi(t) xi(t);
        assume([in(t, 'real'), in(phi(t), 'real'), in(xi(t), 'real')]);    
        O1 = tensor([xi(t),0,0]); 
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);
        S1 = frame(O1,B1); 

        S1_ = S1.subs({t,phi,xi},{2,2*t,3*t});
    end
end
 
function test_plot(~) % tensor plotting  
    import anakin.*
    S = frame;
    
    f = figure;
    S.plot; 
    close(f);
end











