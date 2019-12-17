classdef STsequence  
    % Diffusion gradient : Stejskal-Tanner Pulsed-Gradient Spin Echo 
    properties           % Sequence Parameters
        Delta;           % unit = ms, separation of the two pulses 
        delta;           % unit = ms, duration of the gradient
        Ttot;
        strength;        % typical value ~ e-8T/mu_m 
        direction;       % list of gradient directions, colums are unit vector
        bvalue;          % Stejskal-Tanner formula (without gamma)
    end
    
    methods
        function dgradient = STsequence(Delta,delta,strength,direction)
            dgradient.Delta     = Delta;
            dgradient.delta     = delta;
            dgradient.Ttot      = Delta + delta;
            dgradient.strength  = strength;
            dgradient.direction = direction;
            dgradient.bvalue    = (strength^2)*(delta^2)*(Delta-delta/3);
        end
        
        function [seq_discrete, dt, ds] = time_discretize(dgrad, N_time)
            % discretized time profile, time step, rms jump size
            load('PhysicalConstants.mat');
            dt = dgrad.Ttot/N_time;
            ds = sqrt(2*D*dt);
            seq_discrete = zeros(1, N_time);
            N_d = round(dgrad.delta/dt);
            N_D = round(dgrad.Delta/dt);
            seq_discrete(1 : N_d) = 1;
            seq_discrete(N_D:N_time) = -1;
        end
        
        function DWsignal = CalcSignal(Ph,dgrad,const)
            % calculate diffusion-weighted signal
            dephase = const.gamma*dgrad.strength*dgrad.direction'*Ph;
            DWsignal = real(mean(exp(1i * dephase)),2);
        end
    end
end

