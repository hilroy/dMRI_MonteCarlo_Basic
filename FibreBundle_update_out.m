function [r_new, ph_new] = FibreBundle_update_out(r,ph,dr,d,R,isgradient)
   % case: current position r is outside cyinder with radius R
   % update spin position and dephase
   % d = halved side length    
   % index = (i, j), current lattice index;  i, j are integers
   % initial index is (0,0)
   r_input = r;
   r_prime = r + dr;               % naive position update
   index = floor(((r/d) + 1)/2);   % determine current lattice location
   index_prime = floor(((r_prime/d)+1)/2);
   center = 2 * d * index;
   r_prime = r_prime - center;     % centered coordinate
   
   if index == index_prime
       % naive position update is inside the current lattice, check for 
       % collision with the outer cylinder (multiple collision impossible)
       if norm(r_prime) > R
           % no collision occurs
           r_new = r_prime + center;
       else
           r = r - center;
           coeff = [dot(dr,dr) 2*dot(r,dr) dot(r,r) - R^2];
           t = min(roots(coeff));
           r = r + t * dr;
           normal = r / R;
           Refl = eye(2) - 2 *( normal * normal');
           dr = Refl * ((1 - t) * dr);
           r_new = r + dr + center; 
       end
   else
       % naive update spills to neighbouring lattice
       % update the lattice index
       r_new = r_prime + center;
   end
   ph_new = ph + (r_input + r_new)* isgradient/2;
end

