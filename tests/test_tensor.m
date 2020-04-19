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
function tests = test_tensor
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator 
    import anakin.* 
    B1m = [0,-1,0;1,0,0;0,0,1];
    B1 = basis(B1m);
    E = -3; % scalar
    a = [1;2;0]; % vector
    T = [1,2,3;4,5,6;7,8,9];
        
    % general
    assert(tensor == tensor(0)); % null tensor
    assert(tensor(tensor) == tensor); % itself
    assert(tensor(tensor,B1) == tensor); % basis input basis
    
    % scalar
    assert(all(tensor(E).components == E)); % scalar definition
    
    % vector
    assert(all(tensor(a).components == a)); % components in canonical tensor basis, column array form
    assert(all(tensor(a').components == a)); % see if it rotates a row vector
    assert(tensor(a,B1) == tensor(B1m*a)); % components in another basis, column array form
    assert(tensor(a',B1) == tensor(B1m*a)); % components in another basis, row array form
    
    % tensor
    assert(all(all(tensor(T).components == T))); % components in canonical tensor basis
    assert(tensor(T,B1) == tensor(B1m*T*inv(B1m))); % components in another basis
end

function test_creator_sym(~) % Call creator with sym arguments
    if license('test','symbolic_toolbox') 
        import anakin.*
        syms t theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        
        B1m = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];        
        B1 = basis(B1m);
        E = -2*theta; % scalar
        a = [cos(theta);sin(theta);0]; % vector
        T = [cos(theta),theta^2,3;4,phi-theta,6;7,t^2,9];
        
        % scalar
        assert(isAlways(tensor(E).components == E)); % scalar definition
        
        % vector
        assert(all(isAlways(tensor(a).components == formula(a)))); % components in canonical tensor basis, column array form
        assert(all(isAlways(tensor(a').components == formula(a)))); % see if it rotates a row vector
        assert(tensor(a,B1) == tensor(B1m*a)); % components in another basis, column array form
        assert(tensor(a',B1) == tensor(B1m*a)); % components in another basis, row array form
        
        % tensor
        assert(all(all(isAlways(tensor(T).components == formula(T))))); % components in canonical tensor basis
        assert(tensor(T,B1) == tensor(B1m*T*inv(B1m))); % components in another basis
    end
end

function test_overloads(~) % overloads
    import anakin.*
    
    a = tensor([1.1,2.2,3.3]); 
    b = tensor([4.4,5.5,6.6]);     
    
    assert(a==a);
    assert(a~=b);
    assert(a + b == tensor([5.5;7.7;9.9]));
    assert(a - b == tensor([-3.3;-3.3;-3.3]));
    assert(+a == a);
    assert(-a == tensor([-1.1,-2.2,-3.3]));
    assert(2 * a == tensor([2.2;4.4;6.6]));
    assert(a * 2 == tensor([2.2;4.4;6.6])); 
    assert(2 \ a == tensor([0.55;1.1;1.65]));
    assert(a / 2 == tensor([0.55;1.1;1.65])); 
    assert(dot(a,b) == 38.72);
    assert(norm(a) == sqrt(16.94));
    assert(cross(a,b) == tensor([-3.63,7.26,-3.63]));
    
    assert(unitvector(a) == tensor([1.1,2.2,3.3]/sqrt(16.94)));
    assert(magnitude(a) == sqrt(16.94));
    assert(angle(a,b) == acos(38.72/(sqrt(16.94)*sqrt(93.17))));
    assert(angle(a,b,cross(b,a)) == -acos(38.72/(sqrt(16.94)*sqrt(93.17))));
     
    if license('test','symbolic_toolbox') 
        syms t x1(t) y1(t) z1(t) x2(t) y2(t) z2(t);
        assume([in(t, 'real'), ...
                in(x1(t), 'real'), in(y1(t), 'real'), in(z1(t), 'real'),...
                in(x2(t), 'real'), in(y2(t), 'real'), in(z2(t), 'real')]);
        
        a = tensor([x1,y1,z1]); 
        b = tensor([x2,y2,z2]);     
    
        assert(a==a);
        assert(a~=b);
        assert(a + b == tensor([x1+x2;y1+y2;z1+z2]));
        assert(a - b == tensor([x1-x2;y1-y2;z1-z2]));
        assert(+a == a);
        assert(-a == tensor([-x1,-y1,-z1]));
        assert(2 * a == tensor(2*[x1,y1,z1]));
        assert(a * 2 == tensor(2*[x1,y1,z1])); 
        assert(2 \ a == tensor([x1,y1,z1]/2));
        assert(a / 2 == tensor([x1,y1,z1]/2)); 
        assert(isAlways(dot(a,b) == x1*x2+y1*y2+z1*z2));
        assert(isAlways(norm(a) == sqrt(x1^2+y1^2+z1^2)));
        assert(cross(a,b) == tensor([y1*z2-z1*y2,z1*x2-x1*z2,x1*y2-y1*x2]));
       
        assert(unitvector(a) == tensor([x1,y1,z1]/sqrt(x1^2+y1^2+z1^2)));
        assert(isAlways(magnitude(a) == sqrt(x1^2+y1^2+z1^2)));
        aaa = dot(a,b)/(norm(a)*norm(b));
        assert(isAlways(angle(a,b) == acos(aaa.components)));
        assert(isAlways(angle(a,b,cross(b,a)) == sign(dot(cross(a,b),cross(b,a)))*acos(aaa.components)));
    
    end
end

function test_components(~) % Call components, x,y,z
    import anakin.*
    B1 = basis([0,1,0],[-1,0,0],[0,0,1]); 
    a = tensor([1,2,3],B1);  
    
    assert(all(a.components == [-2;1;3]));
    assert(a.x(1) == -2);
    assert(a.x(2) == 1);
    assert(a.x(3) == 3);
    assert(all(a.components(B1) == [1;2;3]));
    assert(a.x(1,B1) == 1);
    assert(a.x(2,B1) == 2);
    assert(a.x(3,B1) == 3);
        
    if license('test','symbolic_toolbox') 
        syms t theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);        
        a = tensor([cos(theta),sin(theta),0],B1);  

        assert(isAlways(all(a.components(B1) == [cos(theta);sin(theta);0])));
        assert(isAlways(a.x(1,B1) == cos(theta)));
        assert(isAlways(a.x(2,B1) == sin(theta)));
        assert(isAlways(a.x(3,B1) == 0));
    end
end

function test_isa(~) % isunitary, isparallel, isperpendicular...
    import anakin.*
    
    if license('test','symbolic_toolbox')
        ii = 2;
    else
        ii = 1;
    end
    
    for i = 1:ii
    
        if i == 1
            a = tensor([cos(pi/6),sin(pi/6),0]);
            b = tensor([0,0,2.1]);
            T = tensor([1,1i,2;-1i,3,4;2,4,8]);
        elseif i == 2
            syms t theta(t) phi(t);
            assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);    
            a = tensor([cos(theta),sin(theta),0]); 
            b = tensor([0,0,sin(phi)]);                 
        end
        c = 2*a; 

        assert(isunitvector(a)); 
        assert(~isunitvector(c)); 
        assert(isperpendicular(a,b));
        assert(~isperpendicular(a,c));
        assert(~isparallel(a,b));
        assert(isparallel(a,c)); 
        
        assert(ishermitian(T));
        assert(~ishermitian(T+[0,1,0;0,0,0;0,0,0]));
        assert(~isantihermitian(T));
        assert(isantihermitian(tensor([1i,0,0;0,0,0;0,0,0])));
        assert(isunitary(tensor([0,1,0;1,0,0;0,0,1])));
        assert(~isunitary(tensor([1,0,0;0,2,0;0,0,2])));
        assert(isnormal(tensor([1,1,0;0,1,1;1,0,1])));
        assert(~isnormal(tensor([1,1,1;0,1,1;0,0,1])));
    end
end

function test_subs(~) % Particularize a symbolic tensor
    if license('test','symbolic_toolbox')
        import anakin.*
        syms t;
        syms theta(t) phi(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'),in(xi(t), 'real')]);    
        a = tensor([cos(theta),sin(theta)+xi^2*sin(phi),xi*cos(phi)]); 

        c = a.subs({t,theta,phi,xi},{1,t^2-4,2*t+3,-2}); % call with cell arrays
        c = a.subs([t,theta,phi,xi],[1,t^2-4,2*t+3,-2]); % call with arrays
    end
end

function test_timediff(~) % Time derivative wrt a given basis calling dt
    if license('test','symbolic_toolbox')
        import anakin.*
        syms t;
        syms theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        B0 = basis;
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);    
        a = tensor([cos(theta),sin(theta),0],B1); 

        assert(a.dt(B1) == tensor(diff(theta,1)*[-sin(theta),cos(theta),0],B1)); % time derivative of a wrt B1
        assert(a.dt(B0) == a.dt(B1) + cross(omega(B1,B0),a)); % time derivative of a wrt B0
    end
end

function test_plot(~) % tensor plotting 
    import anakin.*
    a = tensor([1,5,3]); 
    
    f = figure;
    a.plot;
    a.plot('color','b');
    a.plot(tensor([1,1,1]),'color',[0.5,0.5,0.5]);
    close(f);
end







