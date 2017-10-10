classdef ParticleFilterX < matlab.mixin.Copyable
    % ParticleFilterX class
    %
    % Summary of ParticleFilterX:
    % This is a class implementation of a bootstrap Particle Filter.
    %
    % ParticleFilterX Properties:
    %    - config       = structure, with fields:
    %       .k                      = time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .Np                (*1) = number of particles. Default = 1000
    %       .particles         (*2) = particles - (nx x Np [x T]) matrix (3rd dimension T is optional)
    %       .w                      = weights - (1 x Np [x T]) vector. Default = repmat(1/Np,Np) (3rd dimension T is optional)  
    %       .x                      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state (Weighted average avg(w*.particles))
    %       .y                      = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .gen_x0            (*1) = function handle of a procedure that samples from the initial pdf p_x0
    %       .resample_strategy (*)  = resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
    %   
    %   - dyn_model (*)   = Object handle to Dynamic Model Class
    %   - obs_model (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies compulsory properties to instantiate a class object
    %   (*x) Signify valid combinations of optional properties required to instantiate a class object
    %
    % ParticleFilterX Methods:
    %    ParticleFilterX - Constructor method
    %    Predict         - Performs UKF prediction step
    %    Update          - Performs UKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs UKF smoothing on a provided set of estimates
    % 
    % ParticleFilterX Example:

    properties
        config
        dyn_model
        obs_model
    end
    
    methods
        function obj = ParticleFilterX(config, dyn_model, obs_model)
        % ParticleFilterX - Constructor method
        %   
        %   Inputs:
        %       config    |
        %       dyn_model | => Check class help for more details
        %       obs_model |
        %   
        %   Usage:
        %       pf = ParticleFilterX(config, dyn_model, obs_model); 
        %
        %   See also Predict, Update, Iterate, Smooth, resample.
        
            % Validate .Np
            if ~isfield(config,'Np') && isfield(config, 'particles')
                fprintf('[PF] Number of particles missing... Assuming "Np = size(particles,2)"..\n');
                config.Np = size(particles,2);
            elseif ~isfield(config,'Np')
                fprintf('[PF] Number of particles missing... Assuming "Np = 1000"..\n');
                config.Np = 1000;
            end
            
            % Validate .particles
            if (~isfield(config,'particles'))
                fprintf('[PF] Particles not given... Proceeding to generation of initial particles..\n');
                if ~isfield(config,'gen_x0')
                    fprintf('[PF] Function handle to sample from initial pdf not given... Cannot proceed..\n');
                    error('[PF] Please supply either an initial set of particles, or a function handle (gen_0) to allow for generation of initial ones!\n');
                else
                    config.particles = config.gen_x0(config.Np)'; % at time k=1
                    config.w = repmat(1/config.Np, 1, config.Np);
                    fprintf('[PF] Generated %d particles with uniform weights\n',config.Np);
                end
            else
                if size(config.particles,2)~=config.Np
                    error('[PF] Given number of particles (Np) is different that the size of supplied particle list! Aborting..\n');
                end
            end
            
            % Validate .w
            if ~isfield(config,'w');
                fprintf('[PF] Initial set of weights not given... Proceeding to auto initialisation!\n');
                config.w = repmat(1/config.Np, config.Np);
                fprintf('[PF] Uniform weights for %d particles have been created\n', config.Np);
            else
                if (all(config.w ==0))
                    fprintf('[PF] Initial set of weights given as all zeros... Proceeding to auto initialisation!\n');
                    config.w = repmat(1/config.Np, 1, config.Np);
                    fprintf('[PF] Uniform weights for %d particles have been created\n', config.Np);   
                end   
            end
             
            % State estimate based on particles
            obj.config.x(:,1) = sum(repmat(config.w,size(config.particles,1),1).*config.particles,2); 
            
            % Validate .y
            if ~isfield(config,'y')
                fprintf('[PF] No initial observation supplied... Assuming "y = 0"\n');
                config.y = 0;
            end
            
            % Validate .resample_strategy
            if ~isfield(config,'resampling_strategy')
                fprintf('[PF] Resampling strategy not given... Assuming "resampling_strategy = systematic_resampling"..\n');
                config.resampling_strategy = 'systematic_resampling';
            end
            
            % Validate .k
            if (~isfield(config,'k')|| config.k<1)
                fprinf('[PF] Iterator (k) was not initialised properly... Setting "k = 1"..');
                config.k = 1;
            end
            
            obj.config = config;
            
            % Validate dyn_model
            if ~isobject(dyn_model); error('[PF] No dynamic model object provided!'); end
            obj.dyn_model = dyn_model;
            
            % Validate obs_model
            if ~isobject(obs_model); error('[PF] No observation model object provided!'); end
            obj.obs_model = obs_model;
      
        end
        
        function Predict(obj)
        % Predict - Performs bootstrap PF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.config.k = 1; % 1 sec)
        %       pf.Predict();
        %
        %   See also ParticleFilterX, Update, Iterate, Smooth, resample.
            
            % Propagate particles through the dynamic model
            obj.config.particles = obj.dyn_model.sys(obj.config.k, obj.config.particles, obj.dyn_model.sys_noise(obj.config.k, obj.config.Np)); % Simply propagate all particles
            
        end
        
        function Update(obj)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "obj.config.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.config.y = y_new; % y_new is the new measurement)
        %       pf.Update(); 
        %
        %   See also ParticleFilterX, Predict, Iterate, Smooth, resample.
        
            if(size(obj.config.y,2)>1)
                error('[PF] More than one measurement have been provided for update. Use ParticleFilterX.UpdateMulti() function instead!');
            elseif size(obj.config.y,2)==0
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
            end
                             
            % Update particle weights
            wk = obj.config.w .*obj.obs_model.eval_likelihood(obj.config.k, obj.config.y , obj.config.particles);
            % Normalize weight vector
            obj.config.w = wk./sum(wk,2);
                        
            % Resampling
            [obj.config.particles, obj.config.w] = obj.resample(obj.config.particles, obj.config.w, obj.config.resampling_strategy);
            
            % Compute estimated state
            obj.config.x(:,1) = sum(repmat(obj.config.w,size(obj.config.particles,1),1).*obj.config.particles,2);             
        end
        
        function UpdateMulti(obj, assocWeights, LikelihoodMatrix)
        % UpdateMulti - Performs bootstrap PF update step, for multiple measurements
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to
        %                       measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %       LikelihoodMatrix: a (Nm x Np) likelihood matrix, where Nm is the number of measurements and Np is the number of particles.
        %
        %   (NOTE: The measurement "obj.config.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.config.y = y_new; % y_new is the new measurement)
        %       pf.Update(); 
        %
        %   See also ParticleFilterX, Predict, Iterate, Smooth, resample.
            ObsNum = size(obj.config.y,2);  
            if(~ObsNum)
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
                return;
            end
            
            if(~exist('assocWeights','var'))
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            if(~exist('LikelihoodMatrix','var') && isfield(obj.config, 'LikelihoodMatrix'))
                LikelihoodMatrix = obj.config.LikelihoodMatrix;
                
            elseif ~exist('LikelihoodMatrix','var')
                LikelihoodMatrix = obj.obs_model.eval_likelihood(obj.config.k, obj.config.y , obj.config.particles);
            end
            expectedLikelihoodMatrix = sum(LikelihoodMatrix.*repmat(assocWeights(2:end)',1,obj.config.Np),1)/sum(assocWeights(2:end));
                
            % Update particle weights
            wk = obj.config.w .* expectedLikelihoodMatrix;
            
            % Normalize weight vector
            obj.config.w = wk./sum(wk,2);
                        
            % Resampling
            [obj.config.particles, obj.config.w] = obj.resample(obj.config.particles, obj.config.w, obj.config.resampling_strategy);
            
            % Compute estimated state
            obj.config.x(:,1) = sum(repmat(obj.config.w,size(obj.config.particles,1),1).*obj.config.particles,2);             
        end
        
        function Iterate(obj)
        % Iterate - Performs a complete bootstrap PF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" and measurement "obj.config.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.config.k = 1; % 1 sec)
        %       (pf.config.y = y_new; % y_new is the new measurement)
        %       pf.Iterate();
        %
        %   See also ParticleFilterX, Predict, Update, Smooth, resample.
           obj.Predict();
           obj.Update();
        end
        
        function smoothed_estimates = Smooth(obj, filtered_estimates)
        % Smooth - Performs FBS smoothing on a provided set of estimates
        %           (Based on [1])
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of obj.config after each iteration
        %   
        %   Outputs:
        %       smoothed_estimates: a copy of the input (1 x N) cell array filtered_estimates, where the .x and .P fields have been replaced with the smoothed estimates   
        %
        %   (Virtual inputs at each iteration)        
        %           -> filtered_estimates{k}.x          : Filtered state mean estimate at timestep k
        %           -> filtered_estimates{k}.P          : Filtered state covariance estimate at each timestep
        %           -> filtered_estimates{k+1}.x_pred   : Predicted state at timestep k+1
        %           -> filtered_estimates{k+1}.P_pred   : Predicted covariance at timestep k+1
        %           -> smoothed_estimates{k+1}.x        : Smoothed state mean estimate at timestep k+1
        %           -> smoothed_estimates{k+1}.P        : Smoothed state covariance estimate at timestep k+1 
        %       where, smoothed_estimates{N} = filtered_estimates{N} on initialisation
        %
        %   (NOTE: The filtered_estimates array can be accumulated by running "filtered_estimates{k} = ukf.config" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       ukf.Smooth(filtered_estimates);
        %
        %   [1] Mike Klaas, Mark Briers, Nando de Freitas, Arnaud Doucet, Simon Maskell, and Dustin Lang. 2006. Fast particle smoothing: if I had a million particles. In Proceedings of the 23rd international conference on Machine learning (ICML '06). ACM, New York, NY, USA, 481-488.
        %
        %   See also ParticleFilterX, Predict, Update, Iterate.
        
            % Allocate memory
            N                           = length(filtered_estimates);
            smoothed_estimates          = cell(1,N);
            smoothed_estimates{N}       = filtered_estimates{N}; 
            
            % Perform Rauch–Tung–Striebel Backward Recursion
            for k = N-1:-1:1
                summ = zeros(obj.config.Np,obj.config.Np);
                for i = 1:obj.config.Np
                    summ(i,:) = obj.dyn_model.eval(filtered_estimates{k}.k, filtered_estimates{k+1}.particles, filtered_estimates{k}.particles(:,i))';    
                    smoothed_estimates{k}.w(1,i) = filtered_estimates{k}.w(1,i) * sum(smoothed_estimates{k+1}.w .* summ(i,:),2);%/sum(filtered_estimates{k}.w(1,i).*p',2);
                end
                for j = 1:obj.config.Np
                    smoothed_estimates{k}.w(1,i) = smoothed_estimates{k}.w(1,i)/sum(filtered_estimates{k}.w(1,:).*summ(i,:),2);
                end
                %smoothed_estimates{k}.w = smoothed_estimates{k}.w./summ;
                smoothed_estimates{k}.w = smoothed_estimates{k}.w./sum(smoothed_estimates{k}.w);
                smoothed_estimates{k}.x = sum(repmat(smoothed_estimates{k}.w,size(filtered_estimates{k}.particles,1),1).*filtered_estimates{k}.particles,2);
            end
        end        
        
        function [xk, wk, idx] = resample(obj, xk, wk, resampling_strategy)
        % resample - Performs particle resampling
        %   
        %   Inputs:
        %       xk : Initial particle set (nx x Np) matrix, where nx is the dimensionality of the state and Np is the number of particles
        %       wk : Initial weights (1 x Np) matrix
        %       resampling_strategy : resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
        %   
        %   Outputs:   
        %       xk : Resampled particle set (nx x Np) matrix, where nx is the dimensionality of the state and Np is the number of particles
        %       wk : Resampled (Uniform) weights (1 x Np) matrix
        %       idx : Indices of resampled particles (Map to old particles)
        %
        %   Usage:
        %       [xk, wk] = resample(obj, xk, wk, resampling_strategy)
        %       [xk, wk, idx] = resample(obj, xk, wk, resampling_strategy)
        %
        %   See also ParticleFilterX, Predict, Update, Smooth, Iterate.    

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
                  edges = min([0 cumsum(wk)],1); % protect against accumulated round-off
                  edges(end) = 1;                 % get the upper edge exact
                  u1 = rand/Np;
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, ~, idx] = histcounts(u1:1/Np:1, edges);
               otherwise
                  error('Resampling strategy not implemented\n')
            end
            xk = xk(:,idx);                    % extract new particles
            wk = repmat(1/Np, 1, Np);          % now all particles have the same weight
        end
        
        
    end
end