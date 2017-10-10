classdef UParticleFilterX < ParticleFilterX
    % UParticleFilterX class
    %
    % Summary of UParticleFilterX:
    % This is a class implementation of an Unscented Particle Filter.
    %
    % UParticleFilterX Properties:
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
    % UParticleFilterX Methods:
    %    ParticleFilterX - Constructor method
    %    Predict         - Performs UKF prediction step
    %    Update          - Performs UKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs UKF smoothing on a provided set of estimates
    % 
    % UParticleFilterX Example:

    properties
        ukf   % Instance of a UKalmanFilter, used to generate optimal proposal
    end
    
    methods
        function obj = UParticleFilterX(config, dyn_model, obs_model)
        % ParticleFilterX - Constructor method
        %   
        %   Inputs:
        %       config    |
        %       dyn_model | => Check class help for more details
        %       obs_model |
        %   
        %   Usage:
        %       upf = UParticleFilterX(config, dyn_model, obs_model); 
        %
        %   See also Predict, Update, Iterate, Smooth, resample.
        
           % Call SuperClass method
           obj@ParticleFilterX(config, dyn_model, obs_model);
           
           % Instantiate UKF 
           obj.ukf = UKalmanFilterX(config, dyn_model, obs_model);
        end
        
        function Predict(obj)
        % Predict - Performs UPF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       Usage:
        %       (upf.config.k = 1; % 1 sec)
        %       (upf.config.y = y_new; % y_new is the new measurement)
        %       upf.Iterate();
        %
        %   See also ParticleFilterX, Update, Iterate, Smooth, resample.
        
            % Update UKF time and measurement
            obj.ukf.config.k = obj.config.k;
            obj.ukf.config.y = obj.config.y;
            
            % Compute UKF prior mean and covariance
            obj.ukf.config.x = sum(repmat(obj.config.w,size(obj.config.particles,1),1).*obj.config.particles,2);
            obj.ukf.config.P = weightedcov(obj.config.particles',obj.config.w');
            
            % Iterate UKF to obtain Optimal Proposal
            obj.ukf.Iterate();
            
            % Sample from UKF proposal
            obj.config.particles = mvnrnd(repmat(obj.ukf.config.x,1,obj.config.Np)', obj.ukf.config.P)'; % Sample from optimal proposal
        end
        
        function Update(obj)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" and measurement "obj.config.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (upf.config.y = y_new; % y_new is the new measurement)
        %       upf.Update(); 
        %
        %   See also UParticleFilterX, Predict, Iterate, Smooth, resample.
            
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
        %       (upf.config.k = 1; % 1 sec)
        %       (upf.config.y = y_new; % y_new is the new measurement)
        %       upf.Iterate();
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, resample.
         
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
        %   See also UParticleFilterX, Predict, Update, Iterate.
            
            % [=====================]
            % [!!! NOT DEVELOPED !!!] 
            % [=====================]
            error("[UPF] UParticleFilterX Smooth function has not been developed yet!");
            
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
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Call SuperClass method
            [xk, wk, idx] = resample@ParticleFilterX(obj, xk, wk, resampling_strategy);
        end
        
        
    end
end