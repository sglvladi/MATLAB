classdef KalmanFilter_new
    properties
        s
        x_smooth
    end
    methods
        function obj = KalmanFilter_new(prop)
%             if ~isfield(prop,'x'); error('Initial state mean missing'); end
%             if ~isfield(prop,'P'); error('Initial state covariance missing'); end
            if ~isfield(prop,'z'); prop.z = []; end
            if ~isfield(prop,'u'); prop.u=0; end
            if ~isfield(prop,'sys'); prop.sys=eye(width(x)); end
            if ~isfield(prop,'B'); prop.B=0; end
            if ~isfield(prop,'Q'); prop.Q=zeros(width(x)); end
            if ~isfield(prop,'R'); error('Observation covariance missing'); end
            if ~isfield(prop,'obs'); prop.obs=eye(width(x)); end
            obj.s = prop;
        end
        
        function s = Predict(obj, s)
            % Prediction for state vector and covariance:
            s.x = s.sys*s.x + s.B*s.u;
            s.P = s.sys*s.P*s.sys' + s.Q;
            
            s.x_pred = s.x;
            s.P_pred = s.P;
            s.F = s.sys;
            s.H = s.obs;
            s.z_pred = s.obs*s.x;
            s.S =s.obs*s.P*s.obs'+s.R;
            
        end
        
        
        function s = Update(obj, s)
            s.x
            s.z
            % Compute Kalman gain factor:
            S_k =s.obs*s.P*s.obs'+s.R;
            K_k = s.P*s.obs'/(S_k);

            % Correction based on observation:
            s.x = s.x + K_k*(s.z-s.obs*s.x);
            s.P = s.P - K_k*s.obs*s.P;
        end
        
        function s = Iterate(obj,s)
            
            % Prediction for state vector and covariance:
            s.x = s.sys*s.x + s.B*s.u;
            s.P = s.sys*s.P*s.sys' + s.Q;
            
            % Compute Kalman gain factor:
            S_k =s.obs*s.P*s.obs'+s.R;
            K_k = s.P*s.obs'/(S_k);

            % Correction based on observation:
            s.x = s.x + K_k*(s.z-s.obs*s.x);
            s.P = s.P - K_k*s.obs*s.P;
            
        end
        
        function obj = Smooth(obj, s)
            r = length(s);
            x_sm = s(r).x;
            P_sm = s(r).P;
            for j=1:r-1
       
                x_min=s(r-j).sys*s(r-j).x;
                P_min=s(r-j).sys*s(r-j).P*s(r-j).sys'+s(r-j).Q;
                G = s(r-j).P*s(r-j).sys'/(P_min);
                x_sm = s(r-j).x + G*(x_sm-x_min);
                obj.x_smooth=[obj.x_smooth, x_sm];
                P_sm=s(r-j).P + G*(s(r-j+1).P-P_sm)*G';
                
            end
        end
            
    end
end