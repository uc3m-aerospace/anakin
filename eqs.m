function eqs = equations(P,B1) % returns vector of equations of motion projected in basis B1
            p = P.p.dt;
            F = anakin.tensor([0;0;0]); % allocate;
            for i=1:length(P.forces)
                F = F + P.forces{i};
            end
            if ~exist('B1','var')
                B1 = anakin.basis;
            elseif isa(B1,'anakin.frame') 
                B1 = B1.basis; % extract basis
            end
            eqs = sym([0;0;0]); % allocate
            for i=1:3
                p_ = p * B1.e(i);
                F_ = F * B1.e(i);
                eqs(i) = (p_.components == F_.components);
            end                
end
        



function eqs = equations(b,B1) % returns vector of equations of motion projected in basis B1
            p = b.p.dt;
            H = b.H.dt;
            F = anakin.tensor([0;0;0]); % allocate;
            M = anakin.tensor([0;0;0]); % allocate;
            for i=1:length(b.forces)
                F = F + b.forces{i};
                M = M + b.torques{i};
            end
            if ~exist('B1','var')
                B1 = anakin.basis;
            elseif isa(B1,'anakin.frame') 
                B1 = B1.basis; % extract basis
            end
            eqs = sym([0;0;0]);
            for i=1:3
                p_ = p * B1.e(i);
                F_ = F * B1.e(i);
                eqs(i) = (p_.components == F_.components);
            end                
            for i=1:3
                H_ = H * B1.e(i);
                M_ = M * B1.e(i);
                eqs(3+i) = (H_.components == M_.components);
            end      
end
        



function test_eqs(~) % Equations
    import anakin.*
    if license('test','symbolic_toolbox') 
        syms t;
        syms theta(t) xi(t);
        assume([in(t, 'real'), in(theta(t), 'real'), in(xi(t), 'real'),...
                in(diff(theta(t),t), 'real'), in(diff(xi(t),t), 'real')]);    
        c = formula([cos(theta);xi^2;xi]);  
        a = tensor(c);
        mass = scalar(3);
        P = particle(mass,a);
        P.forces = {tensor([1,1,xi])};
        
        eqs = P.equations; 
    end
end