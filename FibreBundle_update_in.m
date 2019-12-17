function [r_new, ph_new] = FibreBundle_update_in(r,ph,dr,R,isgradient)
   % case: current position r is inside cyinder with radius R
   % update spin position and dephase
   r_input = r;
   r_prime = r + dr;
   
   if norm(r_prime) < R
       % no reflection
       r_new = r_prime;
   else 
       % to find the collision point, need to solve a quadratic equation 
       coeff = [dot(dr,dr) 2*dot(r,dr) dot(r,r) - R^2]; % equation (*)
       % take the positive root
       t = roots(coeff);
       t = t (t > 0);
       % overwrite: r = collision point, dr = remaining step (reflected) 
       r = r + t * dr;
       normal = - r / R;
       Refl = eye(2) - 2 * (normal * normal');
       dr = Refl * ((1 - t) * dr);
       r_prime = r + dr;
       while norm(r_prime) > R   % multiple collision confirmed
           % for subsequent collisions, (*) reduced to a linear equation
           t = - 2*dot(r, dr) / dot(dr,dr);
           % repeat previous steps
           r = r + t * dr;
           normal = - r / R;
           Refl = eye(2) - 2 * (normal * normal');
           dr = Refl * ((1 - t) * dr);
           r_prime = r + dr;
       end
       r_new = r_prime;
   end
   ph_new = ph + (r_input + r_new) * isgradient/2; % phase update : trapezoidal rule
end

