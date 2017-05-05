classdef KalmanFilter
    properties
        x
        x_post
        P
        P_post
        z 
        u 
        A 
        B 
        Q 
        R
        H 
    end
    methods
        function obj = KalmanFilter(prop)
            if ~isfield(prop,'x'); error('Initial state mean missing'); end
            if ~isfield(prop,'P'); error('Initial state covariance missing'); end
            if ~isfield(prop,'z'); prop.z = []; end
            if ~isfield(prop,'u'); prop.u=0; end
            if ~isfield(prop,'A'); prop.A=eye(width(x)); end
            if ~isfield(prop,'B'); prop.B=0; end
            if ~isfield(prop,'Q'); prop.Q=zeros(width(x)); end
            if ~isfield(prop,'R'); error('Observation covariance missing'); end
            if ~isfield(prop,'H'); prop.H=eye(width(x)); end
            obj.x = prop.x;
            obj.x_post = obj.x;
            obj.P = prop.P;
            obj.P_post = obj.P;
            obj.z = prop.z;
            obj.u = prop.u;
            obj.A = prop.A;
            obj.B = prop.B;
            obj.H = prop.H;
            obj.Q = prop.Q;
            obj.R = prop.R;
        end
        
        function obj = Predict(obj, A, Q)
            obj.A = A;
            obj.Q = Q;
            obj.x = obj.A*obj.x + obj.B*obj.u;
            obj.P = obj.A*obj.P*obj.A' + obj.Q;
        end
        
        
        function obj = Update(obj, z_in)
            obj.z = z_in;
            S_k = obj.H*obj.P*obj.H'+obj.R;
            K_k = obj.P*obj.H'/(S_k);

            % Correction based on observation:
            obj.x = obj.x + K_k*(obj.z-obj.H*obj.x);
            obj.x_post = [obj.x_post, obj.x];
            obj.P = obj.P - K_k*obj.H*obj.P;
        end
        
        function obj = Iterate(obj, A, Q, z_in)
            obj.A = A;
            obj.Q = Q;
            obj.x = obj.A*obj.x + obj.B*obj.u;
            obj.P = obj.A*obj.P*obj.A' + obj.Q;
            
            obj.z = z_in;
            S_k = obj.H*obj.P*obj.H'+obj.R;
            K_k = obj.P*obj.H'/(S_k);

            % Correction based on observation:
            obj.x = obj.x + K_k*(obj.z-obj.H*obj.x);
            obj.x_post = [obj.x_post, obj.x];
            obj.P = obj.P - K_k*obj.H*obj.P;
            obj.P_post = [obj.P_post, obj.P];
        end
        
        function obj = Smooth(obj)
            for j=1:length(obj.x_post)
                objk
            end
        end
            
    end
end