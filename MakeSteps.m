function steps = MakeSteps(ds, dim, N_time, N_walker, type)
   % ds: rms jump size
   % dim: dimension of the step vectors (1,2,3);
   % type: gaussian, uniform, equisized
   % steps: size = dim * N_time * N_walker
   % first generate directions
   N_steps = N_time * N_walker;
   U = rand(1, N_steps);
   if dim == 1
       dir = U;
       dir(dir < 0.5) = -1;
       dir(dir >= 0.5) = 1;
   elseif dim == 2
       theta = exp(2*1i*pi*U);
       dir = [real(theta); imag(theta)];
   elseif dim == 3
       % mathematically proven in Appendix A 
       Z = U;
       r_xy = sqrt(1 - Z.^2);
       U = rand(1, N_steps);
       XY = r_xy.*exp(2*1i*pi*U);
       dir = [real(XY); imag(XY); Z];
   end
   % then multiply dir by appropriate jump size
   if strcmp(type, 'equisized') == 1
       rho = ds * ones(1, N_steps);
   elseif strcmp(type, 'uniform') == 1
       U = rand(1, N_steps);
       rho = sqrt(dim+2) * ds * nthroot(U,dim);
   elseif strcmp(type, 'gaussian') == 1
       U = chi2rnd(dim,[1 N_steps]);
       rho = ds * sqrt(U);
   end
   steps = bsxfun(@times,dir,rho);
   steps = reshape(steps, [dim N_time N_walker]);
end

