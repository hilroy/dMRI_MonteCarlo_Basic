% define geometry
e_x = [1;0;0];
R_in = 2.4;
R_out = 2.5;
d = 3;
FB = FibreBundle(e_x, R_in, R_out, d);

% generate initial positions
N_walker = 10000;
[r_ini_in, N_in, r_ini_out, N_out] = unif_xy(FB, N_walker);

% define diffusion gradient
Delta = 20;
delta = 10;
strength = 4e-4; % this is a bit higher than typical values
direction = [1;0;0];
dgrad = STsequence(Delta,delta,strength,direction);

% discretization
N_time = 10000;
[dgrad_dis, dt, ds] = time_discretize(dgrad, N_time);
steps = MakeSteps(ds, 3, N_time, N_walker, 'uniform');
steps_in = steps(1:2,:,1:N_in);
steps_out = steps(1:2,:,(N_in + 1):N_walker);
steps_z = steps(3,:,:);

% Random Walk
tic
Ph_in = RW_xy(FB, r_ini_in, steps_in, dt, dgrad_dis, 'in');
Ph_out = RW_xy(FB, r_ini_out, steps_out, dt, dgrad_dis, 'out');
Ph_xy = [Ph_in, Ph_out];

Ph_z = RW_free(steps_z, dt, dgrad_dis);
Ph = [Ph_xy ; Ph_z];

Ph = FibreBundle_rotate(Ph, FB.axis);
toc

% calculate signal
load('PhysicalConstants.mat')
dephase = gamma * dgrad.strength * dgrad.direction' * Ph;
signal = real(mean(exp(1i * dephase)));
histogram(dephase)




