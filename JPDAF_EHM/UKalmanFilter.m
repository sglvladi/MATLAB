classdef UKalmanFilter<KalmanFilter_new
%   UKalmanFilter.m                              Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Class implementation of Unscented Kalman Filter. 
%   
%   Properties: 
%     # structure s:
%           s.x    - Initial state
%           s.P    - Initial state covariance
%           s.z    - Initial measurement
%           s.sys  - Dynamics model equation
%           s.obs  - Observation model equation
%           s.Q    - Process noise covariance
%           s.R    - Observation covariance 
%      
%      # alpha 
%      # kapa
%      # beta
%   
%   Dependencies: KalmanFilter_new.m 
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    properties
        alpha
        kapa
        beta
    end
    
    methods
        function obj = UKalmanFilter(prop, alpha, kapa, beta)
            %prop.sys = @(x) x;
            obj@KalmanFilter_new(prop);
            obj.alpha = alpha;
            obj.kapa = kapa;
            obj.beta = beta;
            %obj.s.sys = @(Dt)(@(x)[x(1)+Dt*x(3);x(2)+Dt*x(4);x(3);x(4)]);
        end
        
        % Augmented prediction (used within JPDA)
        function s = Predict(obj,s)
            % Get size of state and measurements vectors
            n=numel(s.x);                                 %numer of states
            m=2;                                 %numer of measurements                                     %default, tunable
            
            % Compute lambda and scaling factor
            lambda=obj.alpha^2*(n+obj.kapa)-n;          %lambda
            c=n+lambda;                                 %scaling factor
            c=sqrt(c);
            
            % Compute weights for sigma points
            [Wm,Wc] = obj.sigma_weights_fl(n, obj.alpha, obj.beta, lambda, c);
            
            % 'Sample' state sigma points
            X = obj.sigmas(s.x, s.P, c);

            % Pass state sigma points through Unscented Transform
            [s.x_pred,X1,s.P_pred,X2]=obj.ut(s.sys,X,Wm,Wc,n,s.Q);          %unscented transformation of process 
            
            
            % 'Sample' measurement sigma points 
            Z=obj.sigmas(s.x_pred, s.P_pred, c);                      %sigma points around x1
           
            % Pass measurement sigma points through Unscented Transform
            [s.z_pred,Z1,s.S,Z2]=obj.ut(s.obs,X1,Wm,Wc,m,s.R);       %unscented transformation of measurments
            
            %Compute transformed cross-covariance
            s.P12=X2*diag(Wc)*Z2';
        
        end
        
        function s = Iterate2(obj,s)
        
            n = size(s.x,1);
            alpha = 1;
            beta = 0;
            kappa = 3-n;

            lambda = alpha^2 * (n + kappa) - n;        
            WM = zeros(2*n+1,1);
            WC = zeros(2*n+1,1);
            for j=1:2*n+1
                if j==1
                    wm = lambda / (n + lambda);
                    wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
                else
                    wm = 1 / (2 * (n + lambda));
                    wc = wm;
                end
                WM(j) = wm;
                WC(j) = wc;
            end

            %
            % Do the filtering
            %  
            %MM = zeros(size(m,1),length(Y));
            %PP = zeros(size(P,1),size(P,2),length(Y));
            %for k=1:length(Y)

            % Form the sigma points for dynamic model
            A = chol(s.P,'lower');
            SX = [zeros(size(s.x)) A -A];
            SX = sqrt(n + lambda)*SX + repmat(s.x,1,size(SX,2));

            % Propagate through the dynamic model
            L=size(SX,2);
            HX =[];
            for k=1:L
                HX(:,k) = s.sys(SX(:,k));
            end
            
            % Compute the predicted mean and covariance
            x_pred = zeros(size(s.x));
            P_pred = zeros(size(s.P));
            for i=1:size(HX,2)
                x_pred = x_pred + WM(i) * HX(:,i);
            end
            for i=1:size(HX,2)
                P_pred = P_pred + WC(i) * (HX(:,i) - x_pred) * (HX(:,i) - x_pred)';
            end
            P_pred = P_pred + s.Q;

            % Form sigma points for measurement step and
            % propagate throught the measurement model
            A = chol(P_pred)';
            SX = [zeros(size(x_pred)) A -A];
            SX = sqrt(n + lambda)*SX + repmat(x_pred,1,size(SX,2));
            L=size(SX,2);
            HY =[];
            for k=1:L
                HY(:,k) = s.obs(SX(:,k));
            end 

            % Compute the updated mean and covariance
            s.x = zeros(size(HY,1),1);
            S  = zeros(size(HY,1),size(HY,1));
            C  = zeros(size(SX,1),size(HY,1));
            for i=1:size(SX,2)
                s.x = s.x + WM(i) * HY(:,i);
            end
            for i=1:size(SX,2)
                S = S + WC(i) * (HY(:,i) - s.x) * (HY(:,i) - s.x)';
                C = C + WC(i) * (SX(:,i) - x_pred) * (HY(:,i) - s.x)';
            end
            S = S + s.R;

            % Compute the gain and updated mean and covariance  
            K = C/S;
            s.x = s.x + K*(s.z - s.x);
            s.P = P_pred - K*S*K';

%             MM(:,k) = m;
%             PP(:,:,k) = P;
            %end
        end
        
        function s = Iterate(obj,s)
            
            % Get size of state and measurements vectors
            n=numel(s.x);                                 %numer of states
            m=numel(s.z);                                 %numer of measurements                                     %default, tunable
            
            % Compute lambda and scaling factor
            lambda=obj.alpha^2*(n+obj.kapa)-n;          %lambda
            c=n+lambda;                                 %scaling factor
            c=sqrt(c);
            
            % Compute weights for sigma points
            [Wm,Wc] = obj.sigma_weights_fl(n, obj.alpha, obj.beta, lambda, c);
            
            % 'Sample' state sigma points
            X = obj.sigmas(s.x, s.P, c);

            % Pass state sigma points through Unscented Transform
            [x1,X1,P1,X2]=obj.ut(s.sys,X,Wm,Wc,n,s.Q);          %unscented transformation of process 
            %Z=obj.sigmas(x1, P1, c);
            
            % 'Sample' measurement sigma points 
            % Z=obj.sigmas(x1, P1, c);                         %sigma points around x1
            % X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
            
            % Pass measurement sigma points through Unscented Transform
            [z1,Z1,P2,Z2]=obj.ut(s.obs,X1,Wm,Wc,m,s.R);       %unscented transformation of measurments
            
            %Compute transformed cross-covariance
            P12=X2*diag(Wc)*Z2';                        
            
            % Compute Kalman gain
            K=P12*inv(P2);
            
            s.x=x1+K*(s.z-z1);                              %state update
            s.P=P1-K*P2*K';                                %covariance update

           
        end

        function [y,Y,P,Y1]=ut(obj,f,X,Wm,Wc,n,R)
            %Unscented Transformation
            %Input:
            %        f: nonlinear map
            %        X: sigma points
            %       Wm: weights for mean
            %       Wc: weights for covraiance
            %        n: numer of outputs of f
            %        R: additive covariance
            %Output:
            %        y: transformed mean
            %        Y: transformed smapling points
            %        P: transformed covariance
            %       Y1: transformed deviations
            L=size(X,2);
            y=zeros(n,1);
            Y=zeros(n,L);
            for k=1:L
                Y(:,k)=f(X(:,k));       
                y=y+Wm(k)*Y(:,k);       
            end
            Y1=Y-y(:,ones(1,L));
            P=Y1*diag(Wc)*Y1'+R;   
        end

        function X=sigmas(obj,x,P,c)
            %Sigma points around reference point
            %Inputs:
            %       x: reference point
            %       P: covariance
            %       c: coefficient
            %Output:
            %       X: Sigma points
            MSR = chol(c^2*P);
            Y = x(:,ones(1,numel(x)));
            X = [x Y+MSR Y-MSR];
        end

        function X=sigma(obj, x, P, alpha, kapa)
            n = numel(x);
            lambda = alpha^2*(n+kapa)-n;
            c = (n+lambda);
            c = sqrt(c);
            MSR = c*chol(P)';
            Y = x(:,ones(1,numel(x)));
            X = [x Y+MSR Y-MSR];
        end

        function [Wm, Wc]=sigma_weights(obj, n, alpha, beta, lambda, c)
            Wm=[lambda/c^2 0.5/c^2+zeros(1,2*n)];          %weights for means
            sum(Wm);
            %Wm = Wm./sum(Wm)
            Wc=Wm;
            %Wc(1)=Wc(1)+(1-alpha^2+beta);
            %Wc = Wc./sum(Wc)
        end
        
        function [Wm, Wc]=sigma_weights_fl(obj, n, alpha, beta, lambda, c)
            Wm = zeros(1,2*n+1);
            Wc = zeros(1,2*n+1);
            Wm(1,1) = lambda/(n + lambda);
            Wc(1,1) = Wm(1,1) + (1 -alpha^2 + beta);
            Wm(1,2:2*n+1) = 1/(2*(n + lambda));
            Wc(1,2:2*n+1) = Wm(1,2:2*n+1);
        end
    end
    
end