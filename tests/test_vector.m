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
function tests = test_vector
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_creator(~) % Call creator without arguments
    import anakin.*
    B1 = basis([0,1,0],[-1,0,0],[0,0,1]); 
    
    a = vector; % null vector
    a = vector(vector); % itself
    a = vector([1;2;0]); % components in canonical vector basis, column array form
    a = vector([1,2,0]); % components in canonical vector basis, row array form    
    a = vector(1,2,0); % components in canonical vector basis, independently given
    a = vector(vector,B1); % relative vector in another basis
    a = vector([1;2;0],B1); % components in another basis, column array form
    a = vector([1,2,0],B1); % components in another basis, row array form       
    a = vector(1,2,0,B1); % components in another basis, independently given

    if license('test','symbolic_toolbox') 
        syms t theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);

        a = vector([cos(theta);sin(theta);0]); % components in canonical vector basis, column array form
        a = vector([cos(theta),sin(theta),0]); % components in canonical vector basis, row array form
        a = vector([cos(theta);sin(theta);0],B1); % components in another basis, column array form
        a = vector([cos(theta),sin(theta),0],B1); % components in another basis, row array form
        a = vector(cos(theta),sin(theta),0); % components in canonical vector basis, independently given
        a = vector(cos(theta),sin(theta),0,B1); % components in another basis, independently given
    end
end


function test_overloads(~) % components, x,y,z
    import anakin.*
    
    a = vector([1.1,2.2,3.3]); 
    b = vector([4.4,5.5,6.6]);     
    
    assert(a==a);
    assert(a~=b);
    assert(a + b == vector([5.5;7.7;9.9]));
    assert(a - b == vector([-3.3;-3.3;-3.3]));
    assert(+a == a);
    assert(-a == vector([-1.1,-2.2,-3.3]));
    assert(2 * a == vector([2.2;4.4;6.6]));
    assert(a * 2 == vector([2.2;4.4;6.6]));
    assert(2 .* a == vector([2.2;4.4;6.6]));
    assert(a .* 2 == vector([2.2;4.4;6.6]));
    assert(2 \ a == vector([0.55;1.1;1.65]));
    assert(a / 2 == vector([0.55;1.1;1.65]));
    assert(2 .\ a == vector([0.55;1.1;1.65]));
    assert(a ./ 2 == vector([0.55;1.1;1.65]));
    assert(dot(a,b) == 38.72);
    assert(norm(a) == sqrt(16.94));
    assert(cross(a,b) == vector([-3.63,7.26,-3.63]));
    
    assert(dir(a) == vector([1.1,2.2,3.3]/sqrt(16.94)));
    assert(magnitude(a) == sqrt(16.94));
    assert(angle(a,b) == acos(38.72/(sqrt(16.94)*sqrt(93.17))));
    assert(angle(a,b,cross(b,a)) == -acos(38.72/(sqrt(16.94)*sqrt(93.17))));
     
    if license('test','symbolic_toolbox') 
        syms t x1(t) y1(t) z1(t) x2(t) y2(t) z2(t);
        assume([in(t, 'real'), ...
                in(x1(t), 'real'), in(y1(t), 'real'), in(z1(t), 'real'),...
                in(x2(t), 'real'), in(y2(t), 'real'), in(z2(t), 'real')]);
        
        a = vector([x1,y1,z1]); 
        b = vector([x2,y2,z2]);     
    
        assert(a==a);
        assert(a~=b);
        assert(a + b == vector([x1+x2;y1+y2;z1+z2]));
        assert(a - b == vector([x1-x2;y1-y2;z1-z2]));
        assert(+a == a);
        assert(-a == vector([-x1,-y1,-z1]));
        assert(2 * a == vector(2*[x1,y1,z1]));
        assert(a * 2 == vector(2*[x1,y1,z1]));
        assert(2 .* a == vector(2*[x1,y1,z1]));
        assert(a .* 2 == vector(2*[x1,y1,z1]));
        assert(2 \ a == vector([x1,y1,z1]/2));
        assert(a / 2 == vector([x1,y1,z1]/2));
        assert(2 .\ a == vector([x1,y1,z1]/2));
        assert(a ./ 2 == vector([x1,y1,z1]/2));
        assert(isAlways(dot(a,b) == x1*x2+y1*y2+z1*z2));
        assert(isAlways(norm(a) == sqrt(x1^2+y1^2+z1^2)));
        assert(cross(a,b) == vector([y1*z2-z1*y2,z1*x2-x1*z2,x1*y2-y1*x2]));
       
        assert(dir(a) == vector([x1,y1,z1]/sqrt(x1^2+y1^2+z1^2)));
        assert(isAlways(magnitude(a) == sqrt(x1^2+y1^2+z1^2)));
        assert(isAlways(angle(a,b) == acos(dot(a,b)/(norm(a)*norm(b)))));
        assert(isAlways(angle(a,b,cross(b,a)) == sign(dot(cross(a,b),cross(b,a)))*acos(dot(a,b)/(norm(a)*norm(b)))));
    
    end
end

function test_components(~) % Call components, x,y,z
    import anakin.*
    B1 = basis([0,1,0],[-1,0,0],[0,0,1]); 
    a = vector([1,2,3],B1);  
    
    assert(all(a.components(B1) == [1;2;3]));
    assert(a.x(B1) == 1);
    assert(a.y(B1) == 2);
    assert(a.z(B1) == 3);
    assert(all(a.components == [-2;1;3]));
    assert(a.x == -2);
    assert(a.y == 1);
    assert(a.z == 3);
    
    if license('test','symbolic_toolbox') 
        syms t theta(t) phi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);
        B1 = basis([1,0,0],[0,cos(phi),sin(phi)],[0,-sin(phi),cos(phi)]);        
        a = vector([cos(theta),sin(theta),0],B1);  

        assert(isAlways(all(a.components(B1) == [cos(theta);sin(theta);0])));
        assert(isAlways(a.x(B1) == cos(theta)));
        assert(isAlways(a.y(B1) ==  sin(theta)));
        assert(isAlways(a.z(B1) == 0));
    end
end

function test_isa(~) % isunitary, isparallel, isperpendicular
    import anakin.*
    
    if license('test','symbolic_toolbox')
        ii = 2;
    else
        ii = 1;
    end
    
    for i = 1:ii
    
        if i == 1
            a = vector(cos(pi/6),sin(pi/6),0);
            b = vector(0,0,2.1);
        elseif i == 2
            syms t theta(t) phi(t);
            assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real')]);    
            a = vector([cos(theta),sin(theta),0]); 
            b = vector([0,0,sin(phi)]);                 
        end
        c = 2*a; 

        assert(isunitary(a)); 
        assert(~isunitary(c)); 
        assert(isperpendicular(a,b));
        assert(~isperpendicular(a,c));
        assert(~isparallel(a,b));
        assert(isparallel(a,c)); 
    end
end

function test_subs(~) % Particularize a symbolic vector
    if license('test','symbolic_toolbox')
        import anakin.*
        syms t;
        syms theta(t) phi(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(phi(t), 'real'),in(xi(t), 'real')]);    
        a = vector([cos(theta),sin(theta)+xi^2*sin(phi),xi*cos(phi)]); 

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
        a = vector([cos(theta),sin(theta),0],B1); 

        assert(a.dt(B1) == vector(diff(theta,1)*[-sin(theta),cos(theta),0],B1)); % time derivative of a wrt B1
        assert(a.dt(B0) == a.dt(B1) + cross(omega(B1,B0),a)); % time derivative of a wrt B0
    end
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







