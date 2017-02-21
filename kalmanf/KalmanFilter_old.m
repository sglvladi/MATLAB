classdef KalmanFilter_old
    properties
        x 
        P
        z 
        u 
        A 
        B 
        Q 
        R
        H 
    end
    methods
        function obj = KalmanFilter(s)
            if ~isfield(s,'x'); error('Initial state mean missing'); end
            if ~isfield(s,'P'); error('Initial state covariance missing'); end
            if ~isfield(s,'z'); error('Observation vector missing'); end
            if ~isfield(s,'u'); s.u=0; end
            if ~isfield(s,'A'); s.sys=eye(width(x)); end
            if ~isfield(s,'B'); s.B=0; end
            if ~isfield(s,'Q'); s.Q=zeros(width(x)); end
            if ~isfield(s,'R'); error('Observation covariance missing'); end
            if ~isfield(s,'H'); s.obs=eye(width(x)); end
            obj.x = s.x;
            obj.P = s.P;
            obj.z = s.z;
            obj.u = s.u;
            obj.sys = s.sys;
            obj.B = s.B;
            obj.obs = s.obs;
            obj.Q = s.Q;
            obj.R = s.R;
        end
        
        function obj = Predict(varargin)
            if nargin >= 1 && isa(varargin{1}, 'KalmanFilter')
                if varargin{2} ~= 0
                    obj.sys = varargin{2};
                end
                if varargin{3} ~= 0
                    obj.Q = varargin{3};
                end
            elseif nargin<1
                    error('Function requires either 1 (kf) or 3 (kf, A, Q) arguments. 0 given');
            end
            obj.x
            obj.x = obj.sys*obj.x + obj.B*obj.u;
            obj.P = obj.sys*obj.P*obj.sys + obj.Q;
        end
        
        
        function obj = Update(obj, measurement)
            pred_state_mean = obj.trans_func*obj.p_state(1);
            pred_state_cov = obj.trans_func*obj.p_state(2)*transpose(obj.trans_function) + obj.proc_noise(2);
            
            pred_meas_mean = obj.likelihood_func*pred_state_mean;
            pred_meas_cov = obj.likelihood_func*pred_state_cov*transpose(obj.likelihood_funct) + obj.meas_noise(2);
            
            K = pred_meas_cov*transpose(obj.likelihood_func)*(obj.likelihood_func*pred_meas_cov*transpose(obj.likelihood_func) + obj.meas_noise(2));
            
            cur_state_mean = pred_state_mean + K*(measurement-pred_meas_mean);
            [m,n] = size(K*obj.likelihood_func);
            cur_state_covd = (eye(m,n) - K*obj.likelihood_func)*pred_meas_cov; 
            obj.p_state(1
        end
           
        function r = roundOff(obj)
            r = round([obj.Value],2);
        end
        function r = multiplyBy(obj,n)
            r = [obj.Value] * n;
        end
    end
end