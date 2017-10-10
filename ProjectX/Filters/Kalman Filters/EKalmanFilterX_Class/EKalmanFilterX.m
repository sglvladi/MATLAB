classdef EKalmanFilterX<KalmanFilterX
    % EKalmanFilterX class
    %
    % Summary of EKalmanFilterX:
    % This is a class implementation of a first-order Extended Kalman Filter.
    %
    % EKalmanFilterX Properties:
    %    - config       = structurem with fields:
    %       .k          = time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .x (*)      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state
    %       .P (*)      = Estimated state covariance (P_{k|k}) - (nx x nx) matrix 
    %       .x_pred     = Predicted state mean (x_{k|k-1}) - (nx x 1) column vector
    %       .P_pred     = Predicted state mean (P_{k|k-1}) - (nx x nx) matrix
    %       .y          = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .y_pred     = Predicted measurement mean (H*x_{k|k-1}) - (ny x 1) column vector
    %       .S          = Innovation covariance (S_k) - (ny x ny) column vector
    %       .K          = Kalman Gain (K_k) - (ny x ny) column vector
    %       .Fjac       = Transition matrix Jacobian - (nx x nx) matrix (Computed internally)
    %       .Hjac       = Observation matrix Jacobian - (ny x ny) matrix (Computed internally)
    %   
    %   - dyn_model (*)   = Object handle to Dynamic Model Class
    %   - obs_model (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies properties necessary to instantiate a class object
    %
    % EKalmanFilterX Methods:
    %    EKalmanFilterX  - Constructor method
    %    Predict         - Performs EKF prediction step
    %    Update          - Performs EKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs EKF smoothing on a provided set of estimates
    % 
    % EKalmanFilterX Example:
    %     N = 500; % Simulate 500 seconds/iterations
    % 
    %     % Constant Velocity Model
    %     config_cv.dim = 2;
    %     config_cv.q = 0.1;
    %     CVmodel = ConstantVelocityModel(config_cv);
    % 
    %     % Positional Observation Model
    %     config_meas.dim = 2;
    %     config_meas.r = 5;
    %     obs_model = PositionalObsModel(config_meas);
    % 
    %     % Set initial true target state
    %     s = zeros(config_cv.dim*2,N);
    %     s(:,1) = [0; 0; 0.1; 0.3]; % (Position target at 0,0 with velocity of 1 m/s on each axis
    % 
    %     % Initiate Kalman Filter
    %     config_kf.k = 1;  % Use k as 1 sec Dt interval for CV model
    %     config_kf.x = s(:,1); 
    %     config_kf.P = CVmodel.config.Q(1);
    %     kf = KalmanFilterX(config_kf, CVmodel, obs_model);
    % 
    %     % Containers
    %     filtered_estimates = cell(1,N);
    % 
    % 
    %     % START OF SIMULATION
    %     % ===================>
    %     for k = 1:N
    % 
    %         % Generate new state and measurement
    %         if(k~=1)
    %             s(:,k) = CVmodel.propagate_parts(s(:,k-1));
    %         end
    %         y(:,k) = obs_model.sample_from(obs_model.transform_mean(s(:,k)));
    % 
    %         % Iterate Kalman Filter
    %         kf.config.y = y(:,k);
    %         kf.Predict(); 
    %         kf.Update();
    % 
    %         % Store filtered estimates
    %         filtered_estimates{k} = kf.config;
    %     end
    %     % END OF SIMULATION
    %     % ===================>
    % 
    %     % Compute smoothed estimates
    %     smoothed_estimates = kf.Smooth(filtered_estimates);
    % 
    %     % Extract estimates
    %     x_filtered = zeros(config_cv.dim*2,N);
    %     x_smoothed = zeros(config_cv.dim*2,N);
    %     for k = 1:N
    %         x_filtered(:,k) = filtered_estimates{k}.x;
    %         x_smoothed(:,k) = smoothed_estimates{k}.x;
    %     end
    % 
    %     figure
    %     plot(s(1,:),s(2,:),'.-k', x_filtered(1,:), x_filtered(2,:), 'o-b', x_smoothed(1,:), x_smoothed(2,:), 'x-r');
    
    properties
    end
    
    methods
        function obj = EKalmanFilterX(config, dyn_model, obs_model)
        % EKalmanFilterX - Constructor method
        %   
        %   Inputs:
        %       config    |
        %       dyn_model | => Check class help for more details
        %       obs_model |
        %   
        %   Usage:
        %       ekf = EKalmanFilterX(config, dyn_model, obs_model); 
        %
        %   See also Predict, Update, Iterate, Smooth.
        
           % Call SuperClass method
            obj@KalmanFilterX(config, dyn_model, obs_model);
        end
        
        function Predict(obj)
        % Predict - Performs EKF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.config.k = 1; % 1 sec)
        %       ekf.Predict();
        %
        %   See also EKalmanFilterX, Update, Iterate, Smooth.
        
            % Prediction for state vector and covariance:
            [obj.config.x_pred,obj.config.Fjac]=obj.jaccsd(@(x)obj.dyn_model.sys(obj.config.k,x), obj.config.x);    %nonlinear update and linearization at current state
            obj.config.P_pred=obj.config.Fjac*obj.config.P*obj.config.Fjac'+obj.dyn_model.config.Q(obj.config.k);                 %partial update

            % Prediction for measurement vector and covariance
            [obj.config.y_pred,obj.config.Hjac]=obj.jaccsd(@(x) obj.obs_model.obs(obj.config.k, x), obj.config.x_pred);    %nonlinear measurement and linearization
            obj.config.S = obj.config.Hjac*obj.config.P_pred*obj.config.Hjac'+ obj.obs_model.config.R(obj.config.k);
        
        end   
        
        function Update(obj)
        % Update - Performs EKF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "obj.config.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.config.y = y_new; % y_new is the new measurement)
        %       ekf.Update(); 
        %
        %   See also EKalmanFilterX, Predict, Iterate, Smooth.
        
            % Compute Kalman gain
            obj.config.K = obj.config.P_pred*obj.config.Hjac'/obj.config.S;
            
            % Compute filtered estimates
            obj.config.x = obj.config.x_pred + obj.config.K*(obj.config.y - obj.config.y_pred);
            obj.config.P = obj.config.P_pred - obj.config.K*obj.config.S*obj.config.K';
        
        end
        
        function UpdateMulti(obj, assocWeights)
        % UpdateMulti - Performs KF update step, for multiple measurements
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to measurements. Default = [0, ones(1,ObsNum)/ObsNum];
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
            ObsDim = size(obj.config.y,1); 
            
            if(~ObsNum)
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                obj.config.x = obj.config.x_pred;
                obj.config.P = obj.config.P_pred;
                return;
            end
            
            if(~exist('assocWeights','var'))
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            
            % Compute Kalman gain
            innov_err      = obj.config.y - obj.config.y_pred(:,ones(1,ObsNum)); % error (innovation) for each sample
            obj.config.K   = obj.config.P_pred*obj.obs_model.config.h(obj.config.k)'/obj.config.S;  

            % update
            %Pc              = (eye(size(obj.dyn_model.config.f(obj.config.k),1)) - obj.config.K*obj.obs_model.config.h(obj.config.k)*obj.config.P_pred);
            Pc              = obj.config.P_pred - obj.config.K*obj.config.S*obj.config.K';
            tot_innov_err   = innov_err*assocWeights(2:end)';
            Pgag            = obj.config.K*((innov_err.*assocWeights(ones(ObsDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*obj.config.K';
            
            obj.config.x    = obj.config.x_pred + obj.config.K*tot_innov_err;  
            obj.config.P    = assocWeights(1)*obj.config.P_pred + (1-assocWeights(1))*Pc + Pgag;
        end
        
        function Iterate(obj)
        % Iterate - Performs a complete EKF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.config.k" and measurement "obj.config.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.config.k = 1; % 1 sec)
        %       (ekf.config.y = y_new; % y_new is the new measurement)
        %       ekf.Iterate();
        %
        %   See also EKalmanFilterX, Predict, Update, Smooth.
        
           % Call SuperClass method
            Iterate@KalmanFilterX(obj);
        end
        
        function smoothed_estimates = Smooth(obj, filtered_estimates)
        % Smooth - Performs EKF smoothing on a provided set of estimates
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of obj.config after each iteration
        %                            
        %   (NOTE: The filtered_estimates array can be computed by running "filtered_estimates{k} = ekf.config" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       ekf.Smooth(filtered_estimates);
        %
        %   See also EKalmanFilterX, Predict, Update, Iterate.
        
            % Allocate memory
            N                           = length(filtered_estimates);
            smoothed_estimates          = cell(1,N);
            smoothed_estimates{N}       = filtered_estimates{N}; 
            
            % Perform Rauch–Tung–Striebel Backward Recursion
            for k = N-1:-1:1
                smoothed_estimates{k}.C     = filtered_estimates{k}.P * filtered_estimates{k+1}.Fjac' / filtered_estimates{k+1}.P_pred;
                smoothed_estimates{k}.x     = filtered_estimates{k}.x + smoothed_estimates{k}.C * (smoothed_estimates{k+1}.x - filtered_estimates{k+1}.x_pred);
                smoothed_estimates{k}.P     = filtered_estimates{k}.P + smoothed_estimates{k}.C * (smoothed_estimates{k+1}.P - filtered_estimates{k+1}.P_pred) * smoothed_estimates{k}.C';                            
            end
        end
        
        function [z,A]=jaccsd(obj, fun,x)
            % JACCSD Jacobian through complex step differentiation
            % [z J] = jaccsd(f,x)
            % z = f(x)
            % J = f'(x)
            %
            z=fun(x);
            n=numel(x);
            h=n*eps;
            A = imag(fun(repmat(x,1,n)+eye(n)*h*1i))./h;
%             for k=1:n
%                 x1=x;
%                 x1(k)=x1(k)+h*1i;
%                 A(:,k)=imag(fun(x1))/h;
%             end
        end
    end
    
end

