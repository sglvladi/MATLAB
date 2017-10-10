classdef PositionalObsModel < handle
    % PositionalObsModel class
    %
    % Summary of PositionalObsModel
    % This is a class implementation of a time & state invariant Linear-Gaussian Observation Model.
    %
    % PositionalObsModel Properties:
    %    config.dim             = dimensionality (1-3, default 2)
    %    config.r               = Observation noise standard deviation (default 0.1 m)
    %    config.h(~)            = Observation matrix handle (Defined internally, but can be overloaded)
    %    config.R(~)            = Observation noise covariance handle (Defined internally, but can be overloaded)
    %
    % PositionalObsModel Methods:
    %    obs                - State-space to measurement-space state projection function h(x_k,v_k)
    %    obs_cov            - State-space to measurement-space covariance projection function h(P_k,R_k)
    %    obs_noise          - Observation noise sample generator
    %    sample             - Sample from measurement model, given a mean and number of samples
    %    eval_likelihood    - Evaluates the likelihood p(y_k|x_k) = N(y_k; x_k, R) of a set of measurements, given a set of (particle) states  

    properties
        config
    end
    
    methods
        function obj = PositionalObsModel(config)
            % Validate .dim
            if ~isfield(config,'dim')
                fprintf('[POModel] Model dimensionality has not been specified... Applying default setting "dim = 2"..\n');
                config.dim = 2;
            end
            
            % Validate .h
            if ~isfield(config,'h')
                switch(config.dim)
                    case(1)
                        config.h = @(~) [1 0]; 
                        %fprintf('[POModel] Observation matrix missing... Applying default setting "h = [1 0]"..\n');
                    case(2)
                        config.h = @(~) [1 0 0 0;
                                         0 1 0 0]; 
                        %fprintf('[POModel] Observation matrix missing... Applying default setting "h = [1 1 0 0]"..\n');
                    case(3)
                        config.h = @(~) [1 0 0 0 0 0;
                                         0 1 0 0 0 0;
                                         0 0 1 0 0 0]; 
                        %fprintf('[POModel] Observation matrix missing... Applying default setting "h = [1 1 1 0 0 0]"..\n');
                end
            end
            
            % Validate .r
            if ~isfield(config,'r')
                config.r = 0.1;
                fprintf('[POModel] Observation noise standard deviation missing... Applying default setting "r = 0.1"..\n');
            end
            
            % Validate .R
            if ~isfield(config,'R')
                switch(config.dim)
                    case(1)
                        config.R = @(~) config.r^2; 
                                    
                        %fprintf('[POModel] Observation noise covariance missing... Applying default setting "R = config.r^2"..\n');
                    case(2)
                        config.R = @(~) eye(2)*config.r^2;
                        %fprintf('[POModel] Observation noise covariance missing... Applying default setting "R = eye(2)*config.r^2"..\n');
                    case(3)
                        config.R = @(~) eye(3)*config.r^2;
                        %fprintf('[POModel] Observation noise covariance missing... Applying default setting "R = eye(3)*config.r^2"..\n');
                end
            end
            
            obj.config = config;
      
        end
        
        function samples = sample(obj, ~, mu, Ns)
        % sample - Sample from measurement model, given a mean and number of samples  
        %
        %   Inputs:
        %       mu: a (ny x 1) mean vector, where ny is the dimensionality of the measurement
        %       Ns: number of samples (Optional, default = 1)
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       samples: a (ny x Ns) matrix of samples, drawn using the observation model   
        %
        %   Usage:
        %   samples = sample(obj, mu, Ns) generates a (ny x Ns) set of Ns samples using the measurement model, given a mean mu.
        %
        %   See also obs, obs_cov, obs_noise, eval_likelihood.
            if(~exist('Ns', 'var'))
                Ns = 1;
            end
            samples = mvnrnd(repmat(mu', Ns, 1), obj.config.R())';
        end
        
        function y_k = obs(obj, ~, x_k, v_k)
        % obs - State-space to measurement-space state projection function h(x_k,v_k) 
        %
        %   Inputs:
        %       x_k: a (nx x Ns) matrix of Ns state vectors from time k, where nx is the dimensionality of the state
        %       v_k: a (ny x Ns) matrix Ns of measurement noise vectors, where ny is the dimensionality of the measurements (Optional)
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       y_k: a (ny x Ns) matrix of Ns state vectors, which have been projected in the measurement-space   
        %
        %   Usage:
        %   y_k = pos.obs(x_k) projects the (nx x Ns) state vector x_k through the observation model, without the inclusion of noise
        %   y_k = cv.sys(x_k, pos.obs_noise(Ns)) projects the (nx x Ns) particle matrix x_k through the observation model, including randomly generated noise by using the PositionalObsModel.obs_noise function.
        %   y_k = cv.sys(x_k, vk) propagates the (nx x Ns) state vector x_k through the observation model, including an externally computed (ny x Ns) noise matrik vk 
        %
        %   See also obs_cov, obs_noise, sample, eval_likelihood.
            
            % Get dimensions
            Ns = size(x_k,2);
            ny = obj.config.dim;
            
            % 
            if(~exist('v_k','var'))
                v_k = zeros(ny,Ns);
            end
            
            y_k = obj.config.h()*x_k + v_k;
        end
        
        function S_k = obs_cov(obj, ~, P_k, R_k)
        % obs_cov - Measurement (Innovation) covariance function h(P_k,R_k) 
        %
        %   Inputs:
        %       P_k: a (nx x nx) state covariance matrix, where nx is the dimensionality of the state
        %       R_k: a (ny x ny) process noise covariance matrix, where nx is the dimensionality of the observations (Optional)
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       S_k: a (ny x ny) innovation covariance    
        %
        %   Usage:
        %   S_k = pos.obs_cov(P_k) computes the (ny x ny) innovation matrix from the (nx x nx) state covariance P_k, without the inclusion of process noise
        %   S_k = pos.obs_cov(P_k, Q_k) computes the (ny x ny) innovation matrix from the (nx x nx) state covariance P_k, including a process noise covariance R_k 
        %
        %   See also obs, obs_noise, sample, eval_likelihood.
        
            % Get dimensions
            ny = obj.config.dim;
        
            if(~exist('R_k','var'))
                R_k = obj.config.R();
            elseif(strcmp(Q_k, 'false'))
                R_k = zeros(ny);
            end
            S_k = obj.config.h()*P_k*obj.config.h()' + R_k;
        end
        
        function v_k = obs_noise(obj, ~, Ns)
        % obs_noise - Observation noise sample generator 
        %
        %   Inputs:
        %       Ns : The number samples to be generated (Optional, default is 1)
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       v_k: a (ny x Ns) matrix of Ns noise samples, where ny is the dimensionality of the observations   
        %
        %   Usage:
        %   v_k = pos.obs_noise() generates a (ny x 1) noise vector
        %   v_k = pos.obs_noise(Ns) generates a (ny x Ns) noise vector 
        %
        %   See also obs, obs_cov, sample, eval_likelihood.
            if(~exist('Ns', 'var'))
                Ns = 1;
            end
            
            % Get dimensions
            ny = obj.config.dim;
        
            v_k = mvnrnd(zeros(ny,Ns)',obj.config.R())';
        end
        
        function likelihood = eval_likelihood(obj, ~, y_k, x_k)
        % eval_likelihood - Evaluates the likelihood p(y_k|x_k) = N(y_k; x_k, R) of a set of measurements, given a set of (particle) state vectors  
        % 
        %   Inputs:
        %       y_k : a (ny x Nm) matrix of Nm observation vectors, where ny is the dimensionality of the observations
        %       x_k : a (nx x Ns) matrix of Ns state vectors, where nx is the dimensionality of the states
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       likelihood: a (Nm x Ns) matrix of likelihoods p(y_k|x_k)    
        %
        %   Usage:
        %   likelihood = eval_likelihood(obj, ~, y_k, x_k) Evaluates and produces a (Nm x Ns) likelikood matrix.
        %
        %   See also obs, obs_cov, obs_noise, sample.
        
            transformed_x_k = obj.config.h()*x_k;
            likelihood = zeros(size(y_k,2), size(x_k,2));
            if(issymmetric(obj.config.R()) && size(x_k,2)>size(y_k,2))
                % If R is symmetric and the number of state vectors is higher than the number of measurements,
                %   then evaluate N(x_k; y_k, R) = N(y_k; x_k, R) to increase speed
                for i=1:size(y_k,2)
                    likelihood(i,:) = mvnpdf(transformed_x_k', y_k(:,i)', obj.config.R())';
                end
            else
                for i=1:size(transformed_x_k,2)
                     likelihood(:,i) = mvnpdf(y_k', transformed_x_k(:,i)', obj.config.R())';  
                end
             end
                        
        end
    end
end