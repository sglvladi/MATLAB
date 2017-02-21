classdef EKalmanFilter<KalmanFilter_new
    %EKALMANFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = EKalmanFilter(prop)
            %prop.sys = @(x) x;
            obj@KalmanFilter_new(prop);
            %obj.s.sys = @(Dt)(@(x)[x(1)+Dt*x(3);x(2)+Dt*x(4);x(3);x(4)]);
        end
        
        function s = Predict(obj,s)
            % Prediction for state vector and covariance:
            [s.x_pred,s.F]=obj.jaccsd(s.sys,s.x);    %nonlinear update and linearization at current state
            s.P_pred=s.F*s.P*s.F'+s.Q;                 %partial update

            % Prediction for measurement vector and covariance
            [s.z_pred,s.H]=obj.jaccsd(s.obs,s.x_pred);    %nonlinear measurement and linearization
            s.S=s.H*s.P_pred*s.H'+ s.R;
        
        end        
        function s = Iterate(obj,s)
            
            % Prediction for state vector and covariance:
            [s.x_pred,s.F]=obj.jaccsd(s.sys,s.x);    %nonlinear update and linearization at current state
            s.P_pred=s.F*s.P*s.F'+s.Q;                 %partial update

            % Prediction for measurement vector and covariance
            [s.z_pred,s.H]=obj.jaccsd(s.obs,s.x_pred);    %nonlinear measurement and linearization
            s.S=s.H*s.P_pred*s.H'+ s.R;
            
            % Compute Kalman gain factor:
            s.K=s.P_pred*s.H'/s.S;
            
            % Correction based on observation:
            s.x=s.x_pred+s.K*(s.z-s.z_pred);
            s.P=s.P_pred-s.K*s.S*s.K';
        end
        
        function obj = Smooth(obj, s)
            r = length(s);
            x_sm = s(r).x;
            P_sm = s(r).P;
            for j=1:r-1
       
                % Prediction for state vector and covariance:
                [x_min,F]=obj.jaccsd(s(r-j).sys,s(r-j).x);    %nonlinear update and linearization at current state
                P_min=F*s(r-j).P*F'+s(r-j).Q;                 %partial update
                P_min = abs(P_min);
                G = s(r-j).P*F'/(P_min);
                x_sm = s(r-j).x + G*(x_sm-x_min);
                obj.x_smooth=[obj.x_smooth, x_sm];
                P_sm=s(r-j).P + G*(s(r-j+1).P-P_sm)*G';
                
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
            m=numel(z);
            A=zeros(m,n);
            h=n*eps;
            for k=1:n
                x1=x;
                x1(k)=x1(k)+h*1i;
                A(:,k)=imag(fun(x1))/h;
            end
        end
    end
    
end

