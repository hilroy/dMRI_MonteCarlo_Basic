function Ph = RW_free(steps, dt, seq_dis)
   % free diffusion, no need for collision checking
   % steps : dim * N_time * N_walker
   traj = cumsum(steps,2);
   dph = dt * bsxfun(@times,traj,seq_dis);
   Ph = trapz(dph,2);
   Ph = reshape(Ph, [size(steps,1), size(steps,3)]);
end

