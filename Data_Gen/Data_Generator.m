classdef Data_Generator
    properties
        s
    end
    
    methods
        
        function obj = Data_Generator(s)
            if nargin<1
                s.sys_dim = 1;
                s.obs_dim = 1;
                s.sys = @(x) x;
                s.obs = @(x) x;
                s.Q = 1;
                s.R = 1;
                s.trackNum = 1;
                s.iterNum = 1;
            end
            obj.s = s;
        end
        
        function obj = setProperties(obj, sys_dim, obs_dim, sys, Q, obs, R, trackNum, iterNum)
            obj.s.sys_dim = sys_dim;
            obj.s.obs_dim = obs_dim;
            obj.s.sys = sys;
            obj.s.obs = obs;
            obj.s.Q = Q;
            obj.s.R = R;
            obj.s.trackNum = trackNum;
            obj.s.iterNum = iterNum;
            
        end
        
        function [tracks] = Generate(obj, s)
            % Check if s is given
                if nargin < 2
                    s = obj.s;
                end
            
            %Get dimensionality of x_init
            sys_dim = obj.s.sys_dim;
            obs_dim = obj.s.obs_dim;
            
            % Initiate tracks list
            tracks = [];
            
            % Define trackObj structure
            trackObj.x_true = zeros(sys_dim, s.iterNum);
            trackObj.observations = zeros(obs_dim, s.iterNum); 
            
            for j = 1:s.trackNum
                s.x_init = [5*rand(); 5*rand(); -5+(10*rand()); -5+(10*rand())];
                % Initiate outputs
                tracks{j}.trackObj = trackObj;

                % Set first x and y
                tracks{j}.trackObj.x_true(:,1) = s.x_init;
                tracks{j}.trackObj.observations(:,1) = s.obs(tracks{j}.trackObj.x_true(:,1)) + mvnrnd(zeros(obs_dim,1),s.R);

                for i=2:s.iterNum
                    tracks{j}.trackObj.x_true(:,i) = s.sys(tracks{j}.trackObj.x_true(:,i-1), i) + mvnrnd(zeros(sys_dim,1),s.Q);
                    tracks{j}.trackObj.observations(:,i) = s.obs(tracks{j}.trackObj.x_true(:,i)) + mvnrnd(zeros(obs_dim,1),s.R);
                end
            end
        end
    end
    
end