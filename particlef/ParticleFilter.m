classdef ParticleFilter
% =====================================================================================
% Parameters:
% pf: structure with the following fields
%       
%       * Variables
%       -------------------
%       .Np               = number of particles
%       .w                = weights   (Np x T)
%       .particles        = particles (dim x Np x T)
%       .xhk              = estimated state
%       .z                = observation vector at time k (column vector)
%       .k                = iteration number
%       
%       * Handles
%       -------------------
%       .sys              = function handle to process equation
%       .obs              = function handle of the observation likelihood PDF p(y[k] | x[k])
%       .sys_noise        = function handle of a procedure that generates system noise
%       .gen_x0           = function handle of a procedure that samples from the initial pdf p_x0
%       .resample_strategy = resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
%

    properties
        pf
    end
    
    methods
        function obj = ParticleFilter(prop)
            % Validate .Np
            if ~isfield(prop,'Np')
                fprintf('Number of particles missing... Assuming "Np = 1"..\n');
                prop.Np = 1;
            end
            
            % Validate .sys
            if ~isfield(prop,'sys')
                fprintf('Function handle to process equation missing... Assuming "sys = .(x)x"..\n');
                prop.sys = @(x) x;
            end
            
            % Validate .particles
            if (~isfield(prop,'particles')|all(prop.particles ==0))
                fprintf('Particles not given... Proceeding to generation of initial particles..\n');
                if ~isfield(prop,'gen_x0')
                    fprintf('Function handle to sample from initial pdf not given... Cannot proceed..\n');
                    error('Please supply either an initial set of particles, or a function handle (gen_0) to allow for generation of initial ones!\n');
                else
                    for i = 1:prop.Np                          % simulate initial particles
                        prop.particles(:,i,1) = prop.gen_x0(); % at time k=1
                    end   
                    prop.w = repmat(1/prop.Np, prop.Np, 1);
                    fprintf('Generated %d particles with uniform weights\n',i);
                    prop.particles(:,:,1)
                    prop.w
                end
            else
                if size(prop.particles,2)~=prop.Np
                    error('Given number of particles (Np) is different that the size of supplied particle list! Aborting..\n');
                end
            end
            
            % Validate .w
            if ~isfield(prop,'w');
                fprintf('Initial set of weights not given... Proceeding to auto initialisation!\n');
                prop.w = repmat(1/prop.Np, prop.Np, 1);
                fprintf('Uniform weights for %d particles have been created\n', prop.Np);
            else
                if (all(prop.w ==0))
                    fprintf('Initial set of weights given as all zeros... Proceeding to auto initialisation!\n');
                    prop.w = repmat(1/prop.Np, prop.Np, 1);
                    fprintf('Uniform weights for %d particles have been created\n', prop.Np);
                end   
            end
             
            
            % Validate .obs
            if ~isfield(prop,'obs')
                error('Function handle for observation likelihood PDF p(y[k]|x[k]) (obs) is not given! Aborting...\n');
            end
            
            % Validate .z
            if ~isfield(prop,'z')
                fprintf('No initial observation supplied... Assuming "z = 0"\n');
                prop.z = 0;
            end
            
            % Validate .sys_noise
            if ~isfield(prop,'sys_noise')
                error('Function handle to generate system noise (sys_noise) has not been given... Aborting..\n');
            end
            
            % Validate .resample_strategy
            if ~isfield(prop,'resampling_strategy')
                fprintf('Resampling strategy not given... Assuming "resampling_strategy = systematic_resampling"..\n');
                prop.resampling_strategy = 'systematic_resampling';
            end
            
            % Validate .k
            if (~isfield(prop,'k')|| prop.k<1)
                fprinf('Iterator (k) was not initialised properly... Setting "k = 1"..');
                prop.k = 1;
            end
            
            obj.pf = prop;
        end
        
        function pf = Iterate(obj, pf)
            k = pf.k;
            if k == 1
               error('error: k must be an integer greater or equal than 2');
            end

            % Initialize variables
            Np = pf.Np;                              % number of particles
            nx = size(pf.particles,1);               % number of states

            % Separate memory
            xkm1 = pf.particles(:,:); % extract particles from last iteration;
            wkm1 = pf.w(:, k-1);                     % weights of last iteration
            xk   = zeros(size(xkm1));     % = zeros(nx,Np);
            wk   = zeros(size(wkm1));     % = zeros(Np,1);

            for i = 1:Np
               % Predict (Propagate particles)
               xk(:,i) = pf.sys(k, xkm1(:,i), pf.sys_noise());

               % Calculate new weights according to observation likelihood
               wk(i) = wkm1(i) * pf.obs(k, pf.z, xk(:,i));
            end;
            % Normalize weight vector
            wk = wk./sum(wk);
            
            % Calculate effective sample size: eq 48, Ref 1
            Neff = 1/sum(wk.^2);
            
            % Resampling
            resample_percentage = 0.50;
            Nt = resample_percentage*Np;
            %if Neff < Nt
               disp('Resampling ...')
               [xk, wk] = obj.resample(xk, wk, pf.resampling_strategy);
               % {xk, wk} is an approximate discrete representation of p(x_k | y_{1:k})
            %end
            
            % Compute estimated state
            pf.xhk(:,k) = zeros(nx,1);
            for i = 1:Np;
               pf.xhk(:,k) = pf.xhk(:,k) + wk(i)*xk(:,i);
            end
            
            % Store new weights and particles
            pf.w(:,k) = wk;
            pf.particles(:,:,k) = xk;
        end
        
        % Resampling function
        function [xk, wk, idx] = resample(obj, xk, wk, resampling_strategy)

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
                  u1 = rand/Np;
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, ~, idx] = histcounts(u1:1/Np:1, edges);
               otherwise
                  error('Resampling strategy not implemented\n')
            end;
            xk = xk(:,idx);                    % extract new particles
            wk = repmat(1/Np, 1, Np);          % now all particles have the same weight
        end
    end
end