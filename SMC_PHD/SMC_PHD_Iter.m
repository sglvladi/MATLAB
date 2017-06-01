function [par] = SMC_PHD_Iter(par)

    %% Initiate storage
    particles   = par.particles;            % Particles
    w           = par.w;                    % Weights
    Np          = par.Np;                   % Number of particles
    sys         = par.sys;                  % Transition function (Dynamics)
    sys_noise   = par.sys_noise;            % Dynamics noise generator
    obs_model   = par.obs_model;            % Observation model (without noise)
    R           = par.R;                    % Observation covariance
    gen_x0      = par.gen_x0;               % Birth particle sampling function
    Pbirth      = par.Pbirth;               % Probability of birth
    Pdeath      = par.Pdeath;               % Probability of death (1-e_(k|k-1))
    J_k         = par.J_k;                  % Number of birth particles
    x_dim       = size(par.particles,1);    % Dimensionality of the state
    y_dim       = size(par.z,1);            % Dimensionality of measurements
    k           = par.k;                    % Time since last iteration (?t)
    z           = par.z;                    % Measurements vector
    PD          = par.PD;                   % Probability of Detection         
    lambda      = par.lambda;               % Clutter rate per unit volume
    %% Predict
    % Expand number of particles to accomodate for births
    particles = [particles, zeros(4, J_k)]; 
    w = [w, zeros(1, J_k)];
    Np_total = Np + J_k;  
    
    % Generate Np normally predicted particles
    particles(:,1:Np) = sys(k, particles(:,1:Np), sys_noise(particles(:,1:Np))); % Simply propagate all particles
    w(:,1:Np) = (1-Pdeath)* w(:,1:Np);

    % Generate birth particles 
    particles(:,Np+1:end) = gen_x0(J_k)';
    w(:,Np+1:end) = 0.2/J_k;
    
    %% Update
    rhi_i = ones(1,size(z,2)); % Assume all measurements are unused
    C_k = zeros(1,size(z,2));
    trans_particles = obs_model(particles(:,:));
    
    % Compute g(z|x)
    g = zeros(size(trans_particles,2),size(z, 2));
    for i = 1:size(z, 2);
        g(:,i) = mvnpdf(trans_particles', z(:,i)', R);
    end
    
    % Compute C_k(z) Eq. (27) of [1]
    for i = 1:size(z,2)   % for all measurements
        C_k(i) = sum(PD*g(:,i)'.*w,2);
    end
    
    % Calculate pi Eq. (21) of [2]
    pi = sum(PD*g./(ones(Np_total,1)*(lambda+C_k)).*(ones(size(z, 2),1)*w)',1);
    
    % Update weights Eq. (28) of [1]
    w = (1-PD + sum(PD*g./(ones(Np_total,1)*(lambda+C_k)),2))'.*w;
    
    %max(pi)
    find(pi>0.8)
    
    N_k = sum(w,2);
    round(N_k)
    
    [par.particles, par.w] = resample(particles, (w/N_k)', par.resampling_strategy, Np);
    par.w = par.w'*N_k;
end

% Resampling function
function [xk, wk, idx] = resample(xk, wk, resampling_strategy, Np_new)

    Np = length(wk);  % Np = number of particles

    % wk = wk./sum(wk); % normalize weight vector (already done)

    switch resampling_strategy
       case 'multinomial_resampling'
          with_replacement = true;
          idx = randsample(1:Np, Np, with_replacement, wk);
        %{
          THIS IS EQUIVALENT TO:
          edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
          edges(end) = 1;                 % get the upper edge exact
          % this works like the inverse of the empirical distribution and returns
          % the interval where the sample is to be found
          [~, idx] = histc(sort(rand(Np,1)), edges);
        %}
       case 'systematic_resampling'
          % this is performing latin hypercube sampling on wk
          edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
          edges(end) = 1;                 % get the upper edge exact
          u1 = rand/Np_new;
          % this works like the inverse of the empirical distribution and returns
          % the interval where the sample is to be found
          [~, ~, idx] = histcounts(u1:1/Np_new:1, edges);
       otherwise
          error('Resampling strategy not implemented\n')
    end
    xk = xk(:,idx);                    % extract new particles
    wk = repmat(1/Np_new, 1, Np_new)';          % now all particles have the same weight
end