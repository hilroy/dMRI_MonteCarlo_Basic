classdef FibreBundle
    % a fibre bundle modeled as periodically placed infinite cylinders
    % this model is created for basic Monte Carlo simulation
    properties
        axis;       % direction of the axis, a unit vector
        radius_in; 
        radius_out; % fibres have thickness
        hfside;     % half of the lattice side length 
        VolFrac;    % volume fractions: [in between out]
    end
    
    methods
        function FibBund = FibreBundle(direction, R_in, R_out, d)
            FibBund.axis = direction;
            FibBund.radius_in = R_in;
            FibBund.radius_out = R_out;
            FibBund.hfside = d;
            A = [pi*R_in^2 pi*R_out^2 4*d^2];
            A = [A(1) diff(A)];
            FibBund.VolFrac = A/(4*d^2);
        end
        
        function [r_ini_in, N_in, r_ini_out, N_out] = unif_xy(FibBund, N_walker)
            % generate 2d random initial positions in the lattice
            % random walk in the longitudinal direction is unrestricted
            d = FibBund.hfside;
            R_in = FibBund.radius_in;
            R_out = FibBund.radius_out;
            
            r_ini_in = [];
            r_ini_out = [];
            
            n_w = 0;
            while n_w < N_walker
                % acceptance-rejection
                r = d * (2*rand(2,1) -1);
                if norm(r) < R_in
                    r_ini_in = [r_ini_in r];
                    n_w = n_w + 1;
                elseif norm(r) > R_out
                    r_ini_out = [r_ini_out r];
                    n_w = n_w + 1;
                end
            end
            N_in = size(r_ini_in,2);
            N_out = size(r_ini_out,2);
        end
        
        function Ph_xy = RW_xy(FibBund, r_ini, steps, dt, seq_dis, type)
            % random walk simulation in the transverse plane of the fibre
            % bundle, depending on type = "in" or "out" 
            N = size(steps,3);
            T = size(steps,2);
            Ph_xy = zeros(2,N);
            
            if strcmp(type, 'in') == 1
                R = FibBund.radius_in;
                for n = 1 : N                   % to be vectorized...
                    r = r_ini(:,n);
                    ph = Ph_xy(:,n);
                    for t = 1 : T 
                        [r , ph] = FibreBundle_update_in(r,ph,steps(:,t,n),R,seq_dis(t));
                    end
                    Ph_xy(:,n) = dt * ph;
                end
            elseif strcmp(type, 'out') == 1
                R = FibBund.radius_out;
                d = FibBund.hfside;
                for n = 1 : N
                    r = r_ini(:,n);
                    ph = Ph_xy(:,n);
                    for t = 1 : T 
                        [r , ph] = FibreBundle_update_out(r,ph,steps(:,t,n),d,R,seq_dis(t));
                    end
                    Ph_xy(:,n) = dt * ph;
                end
            end
        end
    end
end

