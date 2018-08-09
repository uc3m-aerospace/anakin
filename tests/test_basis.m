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
function tests = test_basis 
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    
    B = basis; % canonical basis
    B = basis(basis); % itself     
    for i = 1:2
        if i == 1
            theta = pi/6;
        elseif license('test','symbolic_toolbox') 
            syms t;
            syms theta(t);
            assume([in(t, 'real'), in(theta(t), 'real')]);
        else 
            break;
        end 
        B = basis(vector([cos(theta),sin(theta),0]),vector([-sin(theta),cos(theta),0]),vector([0;0;1])); % with vectors
        B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]); % with arrays 
        B = basis(cos(theta),sin(theta),0,[-sin(theta),cos(theta),0],[0;0;1]); % xyz,j,k
        B = basis([cos(theta),sin(theta),0],-sin(theta),cos(theta),0,[0;0;1]); % i,xyz,k
        B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],0,0,1); % i,j,xyz
        B = basis([cos(theta),sin(theta),0], -sin(theta),cos(theta),0, 0,0,1); % i,xyz,xyz
        B = basis(cos(theta),sin(theta),0, [-sin(theta),cos(theta),0], 0,0,1); % xyz,j,xyz
        B = basis(cos(theta),sin(theta),0,   -sin(theta),cos(theta),0,[0;0;1]); % xyz,xyz,k
        B = basis(cos(theta),sin(theta),0, -sin(theta),cos(theta),0, 0,0,1); % xyz,xyz,xyz
        
        B = basis(cos(theta),sin(theta),0,vector([-sin(theta),cos(theta),0]),vector([0;0;1])); % xyz,j,k
        B = basis(vector([cos(theta),sin(theta),0]),-sin(theta),cos(theta),0,vector([0;0;1])); % i,xyz,k
        B = basis(vector([cos(theta),sin(theta),0]),vector([-sin(theta),cos(theta),0]),0,0,1); % i,j,xyz
        B = basis(vector([cos(theta),sin(theta),0]), -sin(theta),cos(theta),0, 0,0,1); % i,xyz,xyz
        B = basis(cos(theta),sin(theta),0, vector([-sin(theta),cos(theta),0]), 0,0,1); % xyz,j,xyz
        B = basis(cos(theta),sin(theta),0,   -sin(theta),cos(theta),0,vector([0;0;1])); % xyz,xyz,k 
        
        B = basis([sin(theta/2),0,0,cos(theta/2)]); % quaternions
        B = basis([1,0,0,],theta); % axis, angle
    end    
end

function test_logicals(~) % Call logical methods
    import anakin.*
    for i = 1:2
        if i == 1
            theta = pi/6;
        elseif license('test','symbolic_toolbox') 
            syms t;
            syms theta(t);
            assume([in(t, 'real'), in(theta(t), 'real')]);
        else 
            break;
        end  
        B1 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);

        assert(B1.isorthonormal)
        assert(B1.isrighthanded)
    end
end

function test_overloads(~)
    import anakin.*
    for i = 1:2
        if i == 1
            theta = pi/6;
            phi = pi/5;
        elseif license('test','symbolic_toolbox')
            syms t;
            syms theta(t) phi(t);
            assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        else 
            break;
        end  
        B1 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);
        B2 = basis([cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0;0;1]);
        B3 = basis([cos(theta+phi),sin(theta+phi),0],[-sin(theta+phi),cos(theta+phi),0],[0;0;1]);
        B3b = basis([cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0;0;1],B1); 

        assert(B1~=B2)
        assert(B3==B3b)    
        assert(B1*B2 == B3)
    end
end

function test_representations(~)
    import anakin.*
    for i = 1:2
        if i == 1
            theta = pi/6;
            phi = pi/5;
        elseif license('test','symbolic_toolbox')
            syms t;
            syms theta(t) phi(t);
            assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        else 
            break;
        end  
     
        B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);     

        B.matrix;
        B.i.components;
        B.j.components;
        B.k.components;
        B.rotaxis.components;
        B.rotangle;
        B.quaternions;
    end
end

function test_euler(~) % euler angles
    import anakin.*
    for i = 1:2
        if i == 1
            mypsi = pi/7;
            theta = pi/6;
            phi = pi/5;
        elseif i == 2
            mypsi = 2*pi/3;
            theta = pi/3;
            phi = 2*pi/3; 
        end
        
        B0 = basis; % Example of consecutive rotations with Euler angles
        B1 = B0.rotatez(mypsi);
        B2 = B1.rotatex(theta);
        B3 = B2.rotatez(phi);    
        
        assert(vector(B3.euler) == vector([mypsi,theta,phi])); % they are not vectors, but I call vector to simplify comparision here!
        assert(vector(B3.euler(B0,[3,1,3])) == vector([mypsi,theta,phi])); % they are not vectors, but I call vector to simplify comparision here!
        
        B0 = basis; % Example of consecutive rotations with Euler angles
        B1 = B0.rotatex(mypsi);
        B2 = B1.rotatey(theta);
        B3 = B2.rotatez(phi);    
         
        assert(vector(B3.euler(B0,[1,2,3])) == vector([mypsi,theta,phi])); % they are not vectors, but I call vector to simplify comparision here!
    end
end

function test_omegaalpha(~) % Call angular velocity and acceleration wrt canonical vector basis
    import anakin.*
    if license('test','symbolic_toolbox')
        syms t;
        syms theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        
        B = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);     

        assert(B.omega == vector([ 0;0; diff(theta(t), 1)]));
        assert(B.alpha == vector([ 0;0; diff(theta(t), 2)]));
    end
end

function test_omegaalpha2(~) % Call angular velocity, angular acceleration wrt a second vector basis
    import anakin.*
    if license('test','symbolic_toolbox')
        syms t;
        syms theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);    
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0;-sin(phi);cos(phi)]); 
        B2 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1],B1); % defined wrt B1

        assert(B2.omega(B1) == vector([ 0;0; diff(theta(t), 1)],B1));
        assert(B2.alpha(B1) == vector([ 0;0; diff(theta(t), 2)],B1));    
    end
end

function test_subs(~) % subs
    import anakin.*
    if license('test','symbolic_toolbox')
        syms t;
        syms theta(t);
        assume([in(t, 'real'), in(theta(t), 'real')]);
        B1 = basis([cos(theta),sin(theta),0],[-sin(theta),cos(theta),0],[0;0;1]);

        B1_ = B1.subs({t,theta}, {3,t^2}); 
        assert(B1.subs(theta, pi/6) == basis([cos(pi/6),sin(pi/6),0],[-sin(pi/6),cos(pi/6),0],[0;0;1]));
    end
end

function test_rotated(~) % create rotated bases
    import anakin.*
    for i = 1:2
        if i == 1
            mypsi = pi/7;
            theta = pi/6;
            phi = pi/5;
        elseif license('test','symbolic_toolbox')
            syms t;
            syms mypsi(t) theta(t) phi(t);
            assume([in(t, 'real'), in(mypsi(t), 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        else 
            break;
        end   
        B = basis;

        Bx = B.rotatex(pi/6); % with numeric value
        Bx = B.rotatex(theta); % with symbolic expression
        By = B.rotatey(theta);
        Bz = B.rotatez(theta);

        B0 = basis; % Example of consecutive rotations with Euler angles
        B1 = B0.rotatez(mypsi);
        B2 = B1.rotatex(theta);
        B3 = B2.rotatez(phi);    
    end
end

function test_plot(~) % vector plotting  
    import anakin.*
    B = basis;
    B1 = B.rotatex(pi/6).rotatey(pi/6);
    
    f = figure;
    B.plot;
    B1.plot('color','r');
    B.plot(vector([1,1,1]),'color','b');
    close(f);
end












