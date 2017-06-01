classdef UKalmanFilter<KalmanFilter_new
    %UKALMANFILTER Summary of this class goes here
    %   Detailed explanation goes here
  
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
            %Z=obj.sigmas(x1, P1, c);
            
            % 'Sample' measurement sigma points 
            %Z=obj.sigmas(x1, P1, c);                         %sigma points around x1
            % X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
            
            % Pass measurement sigma points through Unscented Transform
            [s.z_pred,Z1,s.S,Z2]=obj.ut(s.obs,X1,Wm,Wc,m,s.R);       %unscented transformation of measurments
            
            %Compute transformed cross-covariance
            s.P12=X2*diag(Wc)*Z2';
            s.S
        
        end
        
        function s = Iterate(obj,s)
        
            nx = numel(s.x);
            nu = nx;
            ny = numel(s.z);
            nv = ny;

            % Parameters for the UKF
            na = nx + nu + nv;
            alpha = 0.5;
            beta = 2;
            kappa = 0;
            lambda = (alpha^2)*(na + kappa) - na;

            % Weights
            Wm = zeros(1,2*na+1);
            Wc = zeros(1,2*na+1);
            Wm(1,1) = lambda/(na + lambda);
            Wc(1,1) = Wm(1,1) + (1 -alpha^2 + beta);
            Wm(1,2:2*na+1) = 1/(2*(na + lambda));
            Wc(1,2:2*na+1) = Wm(1,2:2*na+1);

            % Allocate memory
            Xa_km1 = zeros(na,2*na+1);
            X_k_km1 = zeros(nx,2*na+1);
            Y_k_km1 = zeros(ny,2*na+1);
            P_k_km1 = zeros(nx,nx);
            Pyy = zeros(ny,ny);
            Pxy = zeros(nx,ny);

            % UKF augmented states
            % if det(P_km1) == 0
            %     P_km1 = eye(nx)*1e-6;
            % end
            %P_km1
            %P_km1 = nearestSPD(P_km1);

            Pa_km1 = blkdiag(s.P, s.Q, s.R);
            Xa_km1_m = [s.x; zeros(nu,1); zeros(nv,1)];



            % Scaling parameters and sigma-points
            [Si,flag] = chol((na + lambda)*Pa_km1, 'lower');
            if flag ~= 0
                SP = nearestSPD((na + lambda)*Pa_km1);
                Si = chol(SP, 'lower');
%                 try
%                     SP = chol(s.P)';
%                 catch
%                     sP = nearestSPD(s.P);
%                     SP = chol(sP)';
%                 end
%                     SQ1 = chol(s.Q(1:nx,1:nu))';
%                    %SQ2 = chol(Q_k(4:6,4:6))';
%                    %SQ3 = chol(Q_k(7:9,7:9))';
%                    % SQ4 = chol(Q_k(10:12,10:12))';
%                    % SQ5 = chol(Q_k(13:15,13:15))';
%                    % SQ6 = chol(Q_k(16:18,16:18))';
%                     SQ = blkdiag(SQ1);%,SQ2,SQ3,SQ4,SQ5,SQ6
%                     SR = chol(s.R)';
%                     Si = blkdiag(SP,SQ,SR);
            end
            Xa_km1_m;
            Si;
            Xa_km1(:,1) = Xa_km1_m;
            Xa_km1(:,2:na+1) = Xa_km1_m*ones(1,na) + Si(:,1:na);
            Xa_km1(:,na+2:2*na+1) = Xa_km1_m*ones(1,na) - Si(:,1:na);
            Xa_km1;

            % Prediction (time update)
            for j = 1:2*na+1
                X_k_km1(1:nx,j) = s.sys(Xa_km1(1:nx,j), Xa_km1(nx+1:nx+nu,j));
                Y_k_km1(1:ny,j) = s.obs(X_k_km1(1:nx,j), Xa_km1(nx+nu+1:na,j));
            end

            x_k_km1 = sum((ones(nx,1)*Wm).*X_k_km1,2);
            y_k_km1 = sum((ones(ny,1)*Wm).*Y_k_km1,2);

            % Measurement update
            for j = 1:2*na+1
                P_k_km1 = P_k_km1 + Wc(1,j)*(X_k_km1(1:nx,j) -x_k_km1)*(X_k_km1(1:nx,j) -x_k_km1)';
                Pyy = Pyy + Wc(1,j)*(Y_k_km1(1:ny,j) -y_k_km1)*(Y_k_km1(1:ny,j) -y_k_km1)';
                Pxy = Pxy + Wc(1,j)*(X_k_km1(1:nx,j) -x_k_km1)*(Y_k_km1(1:ny,j) -y_k_km1)';
            end
            K_k = Pxy/Pyy;
            m_k = x_k_km1 + K_k*(s.z -y_k_km1);
            P_k = P_k_km1 - K_k*Pyy*K_k';

            % Return filtered state and covariance
            s.x = m_k;
            s.P = P_k;
        end
        
        function s = Iterate2(obj,s)
            
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
            Z=obj.sigmas(x1, P1, c);                         %sigma points around x1
            % X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
            
            % Pass measurement sigma points through Unscented Transform
            [z1,Z1,P2,Z2]=obj.ut(s.obs,Z,Wm,Wc,m,s.R);       %unscented transformation of measurments
            
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