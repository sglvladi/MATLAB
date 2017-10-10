classdef ConstantHeadingModelX <  matlab.mixin.Copyable % Handle class with copy functionality
    % ConstantHeadingModelX class
    %
    % Summary of ConstantHeadingModelX
    % This is a class implementation of a 2D nonlinear-Gaussian Constant Heading Dynamic Model [1].
    % The model is described by the following SDEs:
    %
    %   dx = v*cos(h)*dt                      | Position on X-axis (m)
    %   dy = v*sin(h)*dt                      | Position on Y axis (m)
    %   dv = q_vel*dW_t,  W_t~N(0,q_vel^2)    | Absolute speed     (m/s)
    %   dh = q_head*dB_t, B_t~N(0,q_head^2)   | Heading            (rad)
    %
    % [1] P. A. Kountouriotis and S. Maskell, "Maneuvering target tracking using an unbiased nearly constant heading model," 2012 15th International Conference on Information Fusion, Singapore, 2012, pp. 2249-2255.
    %
    % ConstantHeadingModelX Properties:
    %    - config   = structure with fields:
    %    	.q_vel            = Velocity noise standard deviation. Default = 0.01 m/s^2
    %    	.q_head           = Heading noise standard deviation. Default = 0.16 rad/s
    %    	.f(Dt, x_k)      = Process transition function handle f(x_k|x_{k-1}), returns matrix (Defined internally, but can be overloaded)
    %    	.Q(Dt)           = Process noise covariance function handle, returns (nx x nx) matrix (Defined internally, but can be overloaded)
    %
    % ConstantHeadingModelX Methods:
    %    sys        - State vector process transition function f(x_{k-1},w_k) 
    %    sys_cov    - State covariance process transition function f(P_{k-1},Q_k)
    %    sys_noise  - Process noise sample generator 


    properties
        config
    end
    
    methods
        function obj = ConstantHeadingModelX(config)
            
            % Validate .f
            if ~isfield(config,'f')
                config.f = @(Dt, xkm1) [xkm1(1,:)+Dt*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+Dt*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:); xkm1(4,:)];
                fprintf('[CHModel] Process transition function missing... Applying default setting "f = @(Dt, xkm1) [xkm1(1,:)+Dt*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+Dt*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:); xkm1(4,:)]"..\n');
            end
            
            % Validate .q_vel
            if ~isfield(config,'q_vel')
                config.q_vel = 0.01;
                fprintf('[CHModel] Velocity noise standard deviation missing... Applying default setting "q_vel = 0.01"..\n');
            end
            
            % Validate .q_head
            if ~isfield(config,'q_head')
                config.q_head = 0.16;
                fprintf('[CHModel] Heading noise standard deviation missing... Applying default setting "q_head = 0.16"..\n');
            end
            
            % Validate .Q
            if ~isfield(config,'Q')
                % Experimental~~~~~~~~~~~~~~~~~~~~~~>
                % config.Q = @(Dt) [Dt^3/3*q_vel*cos(q_head) 0 Dt^2/2*q_vel 0;
                %                  0 Dt^3/3 0 Dt^2/2;
                %                  Dt^2/2 0 Dt 0;
                %                  0 Dt^2/2 0 Dt];
                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
                config.Q = @(Dt) Dt*diag([config.q_vel^4,config.q_vel^4,config.q_vel^2,config.q_head^2]);
                fprintf('[CHModel] Process noise covariance missing... Applying default setting "Q = (Dt) Dt*diag([0,0,config.q_vel^2,config.q_head^2])"..\n');
            end
            obj.config = config;
      
        end
        
        function x_k = sys(obj, Dt, x_km1, w_k)
        % sys - State vector process transition function f(x_{k-1},w_k) 
        %
        %   Inputs:
        %       Dt : Time interval since last timestep (in seconds)
        %       x_km1: a (nx x Ns) matrix of Ns state vectors from time k-1, where nx is the dimensionality of the state
        %       w_k: a (nx x Ns) matrix Ns of process noise vectors, corresponding to the state vectors x_km1. (Optional)  
        %
        %   Outputs:
        %       x_k: a (nx x Ns) matrix of Ns state vectors, which have been propagated through the dynamic model   
        %
        %   Usage:
        %   xk = cv.sys(Dt, x_km1) propagates the (nx x Ns) state vector x_km1 through the dynamic model for time Dt, without the inclusion of process noise
        %   pk = cv.sys(Dt, p_km1, cv.sys_noise(Ns)) propagates the (nx x Ns) particle matrix p_km1 through the dynamic model for time Dt, including randomly generated process noise by using the ConstantVelocityModel.sys_noise function.
        %   sk = cv.sys(Dt, s_km1, wk) propagates the (nx x Ns) sigma-point matrix s_km1 through the dynamic model for time Dt, including an externally computed (nx x Ns) process noise matrik wk 
        %
        %   See also sys_cov, sys_noise.
        
             % Get dimensions
            Ns = size(x_km1,2);
            nx = size(x_km1,1);
            
            % 
            if(~exist('w_k','var'))
                w_k = zeros(nx,Ns);
            end
            
            x_k = obj.config.f(Dt,x_km1) + w_k;
        end
        
        function P_k = sys_cov(obj, Dt, P_km1, Q_k)
        % sys_cov - State covariance process transition function f(P_{k-1},Q_k) 
        %
        %   Inputs:
        %       Dt : Time interval since last timestep (in seconds)
        %       P_km1: a (nx x nx) state covariance matrix, where nx is the dimensionality of the state
        %       Q_k: a (nx x nx) process noise covariance matrix. (Optional)  
        %
        %   Outputs:
        %       P_k: a (nx x nx) state covariance, which have been propagated through the dynamic model   
        %
        %   Usage:
        %   Pk = cv.sys_cov(Dt, P_km1) propagates the (nx x nx) state covariance P_km1 through the dynamic model for time Dt, without the inclusion of process noise
        %   Pk = cv.sys_cov(Dt, P_km1, Q_k) propagates the (nx x nx) state covariance P_km1 through the dynamic model for time Dt, including a process noise covariance Q_k 
        %
        %   See also sys, sys_noise.
        
            % Get dimensions
            nx = size(P_km1,1);
        
            if(~exist('Q_k','var'))
                Q_k = obj.config.Q(Dt);
            elseif(strcmp(Q_k, 'false'))
                Q_k = zeros(nx);
            end
            P_k = obj.config.f(Dt, P_km1) + Q_k;
        end
        
        function w_k = sys_noise(obj, Dt, Ns)
        % sys_noise - Process noise sample generator 
        %
        %   Inputs:
        %       Dt : Time interval since last timestep (in seconds)
        %       Ns : The number samples to be generated (Optional, default is 1)  
        %
        %   Outputs:
        %       w_k: a (nx x Ns) matrix of Ns process noise samples, where nx is the dimensionality of the state   
        %
        %   Usage:
        %   w_k = cv.sys_noise(Dt) generates a (nx x 1) noise vector
        %   w_k = cv.sys_noise(Dt, Ns) generates a (nx x Ns) noise vector 
        %
        %   See also sys, sys_cov.
            if(~exist('Ns', 'var'))
                Ns = 1;
            end
            
            % Get dimensions
            nx = 4;
        
            w_k = mvnrnd(zeros(nx,Ns)',obj.config.Q(Dt))';
        end
        
    end
end