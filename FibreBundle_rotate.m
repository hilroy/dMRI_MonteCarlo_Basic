function Ph_rot = FibreBundle_rotate(Ph, dir)   
   % where the fibre axis is assumed to be Z
   % dir is the actual fibre direction
   % Ph_rot is the adjusted phase sample, by a 3d rotation
   e_z = [0;0;1];
   k = cross(e_z, dir);                   % axis of rotation
   s = norm(k);                           % sin(angle of rotation)
   if s == 0
       Ph_rot = Ph;
   else
       c = dot(e_z, dir);                 % cos(angle of rotation)
       k_cross = dir * e_z' - e_z * dir'; % matrix notation of cross(k, )
       k_sqd = k * k' /(1 + c);           % "k squared"
       Rot = c*eye(3) + k_cross + k_sqd;  % 3d rotation matrix
       Ph_rot = Rot * Ph;
   end
end

