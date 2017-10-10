classdef EParticleFilterX < ParticleFilterX
    % EParticleFilterX class
    %
    % Summary of EParticleFilterX:
    % This is a class implementation of an Extended Particle Filter.
    %
    % EParticleFilterX Properties:
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
    % EParticleFilterX Methods:
    %    EParticleFilterX - Constructor method
    %    Predict         - Performs EPF prediction step
    %    Update          - Performs EPF update step
    %    Iterate         - Performs a complete EPF iteration (Predict & Update)
    %    Smooth          - Performs EPF smoothing on a provided set of estimates
    % 
    % EParticleFilterX Example:

    properties
        ekf   % Instance of a EKalmanFilter, used to generate optimal proposal
    end
    
    methods
        function obj = EParticleFilterX(config, dyn_model, obs_model)
        % ParticleFilterX - Constructor method
        %   
        %   Inputs:
        %       config    |
        %       dyn_model | => Check class help for more details
        %       obs_model |
        %   
        %   Usage:
        %       epf = EParticleFilterX(config, dyn_model, obs_model); 
        %
        %   See also Predict, Update, Iterate, Smooth, resample.
        
           % Call SuperClass method
           obj@ParticleFilterX(config, dyn_model, obs_model);
           
           % Instantiate EKF
           obj.ekf = EKalmanFilterX(config, dyn_model, obs_model);
        end
        
        function Predict(obj)
        % Predict - Performs EPF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       Usage:
        %       (epf.config.k = 1; % 1 sec)
        %       (epf.config.y = y_new; % y_new is the new measurement)
        %       epf.Iterate();
        %
        %   See also EParticleFilterX, Update, Iterate, Smooth, resample.
        
            % Update EKF time and measurement
            obj.ekf.config.k = obj.config.k;
            obj.ekf.config.y = obj.config.y;
            
            % Compute EKF prior mean and covariance
            obj.ekf.config.x = sum(repmat(obj.config.w,size(obj.config.particles,1),1).*obj.config.particles,2);
            obj.ekf.config.P = weightedcov(obj.config.particles',obj.config.w');
            
            % Iterate EKF to obtain Optimal Proposal
            obj.ekf.Iterate();
            
            % Sample from EKF proposal
            obj.config.particles = mvnrnd(repmat(obj.ekf.config.x,1,obj.config.Np)', obj.ekf.config.P)'; % Sample from optimal proposal
        end
        
        function Update(obj)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" and measurement "obj.config.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (epf.config.y = y_new; % y_new is the new measurement)
        %       epf.Update(); 
        %
        %   See also EParticleFilterX, Predict, Iterate, Smooth, resample.
        
            % Call SuperClass method
            Update@ParticleFilterX(obj);             
        end
        
        function Iterate(obj)
        % Iterate - Performs a complete bootstrap PF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" and measurement "obj.config.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (epf.config.k = 1; % 1 sec)
        %       (epf.config.y = y_new; % y_new is the new measurement)
        %       epf.Iterate();
        %
        %   See also EParticleFilterX, Predict, Update, Smooth, resample.
        
           % Call SuperClass method
           Iterate@ParticleFilterX(obj);
        end
        
        function smoothed_estimates = Smooth(obj, filtered_estimates)
        % Smooth - Performs FBS smoothing on a provided set of estimates
        %   
        %    [=====================]
        %    [!!! NOT DEVELOPED !!!] 
        %    [=====================]
        %
        %   See also EParticleFilterX, Predict, Update, Iterate.
            
            % [=====================]
            % [!!! NOT DEVELOPED !!!] 
            % [=====================]
            error("[EPF] EParticleFilterX Smooth function has not been developed yet!");
            
            % Call SuperClass method 
            % smoothed_estimates = Smooth@ParticleFilterX(obj);
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
        %   See also EParticleFilterX, Predict, Update, Smooth, Iterate.    
        
            % Call SuperClass method
            [xk, wk, idx] = resample@ParticleFilterX(obj, xk, wk, resampling_strategy);
        end
    end
end