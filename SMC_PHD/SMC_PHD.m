classdef SMC_PHD < handle
% =====================================================================================
% Parameters:
% config: structure with the following parameters
%       
%       * Variables
%       -------------------
%       .type               = Type of PHD filter to use. (Currently only 'standard' and 'search' are implemented)
%       .particles          = Particles
%       .w                  = Weights
%       .Np;                = Number of particles
%       .R                  = Observation covariance
%       .Pbirth             = Probability of birth
%       .Pdeath             = Probability of death (1-e_(k|k-1))
%       .PD                 = Probability of Detection  
%       .J_k                = Number of birth particles
%       .z;                 = Measurements vector
%       .lambda;            = Clutter rate per unit volume
%       .resample_strategy  = resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
%       .k                  = Time index or time since last iteration (Dt)
%                             (Depends on definition of sys handle)
%       
%       * Handles
%       -------------------
%       .sys;               = Transition function (Dynamics)
%       .sys_noise;         = Dynamics noise generator
%       .likelihood         = Likelihood pdf p(z|x)
%       .obs_model          = Observation model (without noise)
%       .gen_x0             = Birth particle sampling function
%

    properties
        config
    end
    
    methods
        function obj = SMC_PHD(prop)
            % Validate .Np
            if ~isfield(prop,'Np')
                fprintf('Number of particles missing... Assuming "Np = 1"..\n');
                prop.Np = 1;
            end
            
            % Validate .sys
            if ~isfield(prop,'sys')
                fprintf('Function handle to process equation missing... Assuming "sys = .(x)x"..\n');
                prop.sys = @(x) x;
            end
            
            % Validate .particles
            if (~isfield(prop,'particles')|all(prop.particles ==0))
                fprintf('Particles not given... Proceeding to generation of initial particles..\n');
                if ~isfield(prop,'gen_x0')
                    fprintf('Function handle to sample from initial pdf not given... Cannot proceed..\n');
                    error('Please supply either an initial set of particles, or a function handle (gen_0) to allow for generation of initial ones!\n');
                else
                    %for i = 1:prop.Np                          % simulate initial particles
                        %prop.gen_x0()
                    prop.particles(:,:) = prop.gen_x0(prop.Np)'; % at time k=1
                    %end   
                    prop.w = repmat(1/prop.Np, prop.Np, 1);
                    fprintf('Generated %d particles with uniform weights\n',prop.Np);
                    prop.xhk(:,1) = sum(bsxfun(@times, prop.w(:,1)', prop.particles(:,:)),2);
                end
            else
                if size(prop.particles,2)~=prop.Np
                    error('Given number of particles (Np) is different that the size of supplied particle list! Aborting..\n');
                end
            end
            
            % Validate .w
            if ~isfield(prop,'w');
                fprintf('Initial set of weights not given... Proceeding to auto initialisation!\n');
                prop.w = repmat(1/prop.Np, prop.Np);
                fprintf('Uniform weights for %d particles have been created\n', prop.Np);
            else
                if (all(prop.w ==0))
                    fprintf('Initial set of weights given as all zeros... Proceeding to auto initialisation!\n');
                    prop.w = repmat(1/prop.Np, prop.Np);
                    fprintf('Uniform weights for %d particles have been created\n', prop.Np);
                end   
            end
            
            % Validate .z
            if ~isfield(prop,'z')
                fprintf('No initial observation supplied... Assuming "z = 0"\n');
                prop.z = 0;
            end
            
             % Validate .obs_model
            if ~isfield(prop,'obs_model')
                error('Function handle for observation model (obs_model) (without noise) has not been given... Aborting..\n');
            end
            
             % Validate .likelihood
            if ~isfield(prop,'likelihood')
                error('Function handle for likelihood model p(y|x)(likelihood) has not been given... Aborting..\n');
            end
            
            % Validate .sys_noise
            if ~isfield(prop,'sys_noise')
                error('Function handle to generate system noise (sys_noise) has not been given... Aborting..\n');
            end
            
            % Validate .resample_strategy
            if ~isfield(prop,'resampling_strategy')
                fprintf('Resampling strategy not given... Assuming "resampling_strategy = systematic_resampling"..\n');
                prop.resampling_strategy = 'systematic_resampling';
            end
            
            % Validate .k
            if (~isfield(prop,'k')|| prop.k<1)
                fprinf('Iterator (k) was not initialised properly... Setting "k = 1"..');
                prop.k = 1;
            end
            
            if ~isfield(prop,'Pdeath')
                fprintf('Probability Pdeath missing... Assumming Pdeath = 0.005.\n');
                prop.Pdeath = 0.005;
            end
            
            if ~isfield(prop,'PD')
                fprintf('Probability PD missing... Assumming PD = 0.9.\n');
                prop.PD = 0.9;
            end
            
            if ~isfield(prop,'P_conf')
                fprintf('Probability P_conf missing... Assumming P_conf = 0.9.\n');
                prop.P_conf = 0.9;
            end
            
            if ~isfield(prop,'birth_strategy')
                fprintf('birth_strategy missing... Assumming birth_strategy = "expansion".\n');
                prop.P_conf = 'expansion';
            end
            obj.config = prop;
      
        end
        
        function [config] = Predict(obj)
            switch obj.config.type
                case 'standard'
                    config = obj.Predict_Standard();
                case 'search'
                    config = obj.Predict_Search();
            end
        end
        
        function [config] = Update(obj)
            switch obj.config.type
                case 'standard'
                    config = obj.Update_Standard();
                case 'search'
                    config = obj.Update_Search();
            end
        end
        
        function [config] = Predict_Standard(obj)
            
            % local copy of config
            config = obj.config;
            
            if(strcmp(config.birth_strategy, 'expansion'))
            
                % Expand number of particles to accomodate for births
                config.particles = [config.particles, zeros(4, config.J_k)]; 
                config.w = [config.w, zeros(1, config.J_k)];
                config.Np_total = config.Np + config.J_k;  

                % Generate Np normally predicted particles
                config.particles(:,1:config.Np) = config.sys(config.k, config.particles(:,1:config.Np), config.sys_noise(config.Np)); % Simply propagate all particles
                config.w(:,1:config.Np) = (1-config.Pdeath)* config.w(:,1:config.Np);

                % Generate birth particles 
                config.particles(:,config.Np+1:end) = config.gen_x0(config.J_k)';
                config.w(:,config.Np+1:end) = 0.2/config.J_k;
            elseif(strcmp(config.birth_strategy, 'mixture'))
                a = config.Pbirth;
                b = (1-config.Pdeath);
                Np_n = binornd(config.Np,b/(a+b));

                % Generate normally predicted particles 
                if(Np_n)
                    config.particles(:,1:Np_n) = config.sys(config.k, config.particles(:,1:Np_n), config.sys_noise(config.particles(:,1:Np_n))); % Simply propagate all particles
                end

                % Generate birth particles 
                if(Np_n<config.Np)
                    config.particles(:,Np_n+1:end) = config.gen_x0(config.Np-(Np_n))';
                end

                config.Np_total = config.Np;
            else
                error('Birth strategy "%s" not defined', config.birth_strategy); 
            end
            
            % reassing config
            obj.config = config;
        end
        
        function [config] = Update_Standard(obj)
            %% Update
            
            config = obj.config;
            
             % Tranform particles to measurement space
            trans_particles = config.obs_model(config.particles(:,:)); 
            
            % Compute g(z|x) matrix 
            config.g = zeros(size(trans_particles,2),size(config.z, 2));
            for i = 1:size(config.z, 2)
                config.g(:,i) = config.likelihood(config.k, trans_particles', config.z(:,i)');
            end
            
            % Compute C_k(z) Eq. (27) of [1]  
            C_k = zeros(1,size(config.z,2));
            for i = 1:size(config.z,2)   % for all measurements
                C_k(i) = sum(config.PD*config.g(:,i)'.*config.w,2);
            end
            config.C_k = C_k;
            
            % Update weights Eq. (28) of [1]
            config.w = (1-config.PD + sum(config.PD*config.g./(ones(config.Np_total,1)*(config.lambda+config.C_k)),2))'.*config.w;

            N_k = sum(config.w,2);
            round(N_k)

            [config.particles, config.w] = obj.resample(config.particles, (config.w/N_k)', config.resampling_strategy, config.Np);
            config.w = config.w'*N_k;
            obj.config = config;
        end
        
        function [config] = Predict_Search(obj)
            %% Predict
            
            config = obj.config;
            
            if(strcmp(config.birth_strategy, 'expansion'))
            
                % Expand number of particles to accomodate for births
                config.particles = [config.particles, zeros(4, config.J_k)]; 
                config.w = [config.w, zeros(1, config.J_k)];
                config.Np_total = config.Np + config.J_k;  

                % Generate Np normally predicted particles
                config.particles(:,1:config.Np) = config.sys(config.k, config.particles(:,1:config.Np), config.sys_noise(config.Np)); % Simply propagate all particles
                config.w(:,1:config.Np) = (1-config.Pdeath)* config.w(:,1:config.Np);

                % Generate birth particles 
                config.particles(:,config.Np+1:end) = config.gen_x0(config.J_k)';
                config.w(:,config.Np+1:end) = 0.2/config.J_k;
            elseif(strcmp(config.birth_strategy, 'mixture'))
                a = config.Pbirth;
                b = (1-config.Pdeath);
                Np_n = binornd(config.Np,b/(a+b));

                % Generate normally predicted particles 
                if(Np_n)
                    config.particles(:,1:Np_n) = config.sys(config.k, config.particles(:,1:Np_n), config.sys_noise(config.particles(:,1:Np_n))); % Simply propagate all particles
                end

                % Generate birth particles 
                if(Np_n<config.Np)
                    config.particles(:,Np_n+1:end) = config.gen_x0(config.Np-(Np_n))';
                end
                config.w(:, Np_n+1:end) = 0.2/(config.Np-Np_n);
                config.Np_total = config.Np;
            else
                error('Birth strategy "%s" not defined.. Choose between "expansion" or "mixture" strategies!', config.birth_strategy); 
            end
            
            obj.config = config;
            
        end
        
        function [config] = Update_Search(obj)
            %% Update
            
            config = obj.config;
            
            % Tranform particles to measurement space
            trans_particles = config.obs_model(config.particles(:,:)); 
            
            % Get rhi measurement weights
            config.rhi = config.rhi==1; %ones(1,size(config.z,2)); % Assume all measurements are unused
            
            % Compute g(z|x) matrix (Np x 
            config.g = zeros(size(trans_particles,2),size(config.z, 2));
            for i = 1:size(config.z, 2)
                config.g(:,i) = config.likelihood(config.k, trans_particles', config.z(:,i)');
            end

            % Compute C_k(z) Eq. (27) of [1]  
            C_k = zeros(1,size(config.z,2));
            for i = 1:size(config.z,2)   % for all measurements
                C_k(i) = sum(config.PD*config.rhi(i)*config.g(:,i)'.*config.w,2);
            end
            config.C_k = C_k;

            % Calculate pi Eq. (21) of [2]
            config.pi = zeros(1, size(config.z,2));
            for j = 1:size(config.z,2)
                config.pi(j) = sum((config.PD*config.rhi(j)*config.g(:,j)'/(config.lambda+config.C_k(j))).*config.w,2);
            end
            %config.pi = sum(config.PD*repmat(config.rhi,config.Np_total,1).*config.g./(ones(config.Np_total,1)*(config.lambda+config.C_k)).*(ones(size(config.z, 2),1)*config.w)',1);
            
            % Update weights Eq. (28) of [1]
            w = zeros(size(config.z,2)+1, config.Np_total);
            w(1,:) = (1-config.PD)*config.w;
            for j = 1:size(config.z,2)
                w(j+1,:) = (config.PD*config.rhi(j)*config.g(:,j)'/(config.lambda+config.C_k(j))).*config.w; 
            end
            %w(2:end,:) = (config.PD*repmat(config.pi,config.Np_total,1).*config.g./repmat(config.lambda+config.C_k,config.Np_total,1))'.*repmat(config.w,size(config.pi,2),1);
            
            
            % Select measurements to be used for new tracks (pi>0.9)
            CritMeasurements = find(config.pi>config.P_conf);
            config.NewTracks = [];
            for j = 1:size(CritMeasurements,2)
                MeasInd = CritMeasurements(j); % Index of measurement
                
                % Get particles and weights
                NewTrack.particles = config.particles;
                try
                    NewTrack.w         = w(MeasInd+1,:);
                catch
                    sd =2;
                end
                % Resample particles to ensure they are correctly localised
                % around the measurement
                N_k = sum(NewTrack.w,2);
                [NewTrack.particles, NewTrack.w] = obj.resample(NewTrack.particles, (NewTrack.w/N_k)', config.resampling_strategy, config.Np_conf);
                NewTrack.ExistProb = config.pi(MeasInd);
                config.NewTracks{end+1} = NewTrack; 
            end
            config.pi
            CritMeasurements;
            % Select measurements which are not to be used for new tracks
            NonCritMeasurements = setdiff([1:size(config.z, 2)], CritMeasurements);
            config.w = sum(w([1,NonCritMeasurements+1],:),1);
            %config.w = (1-config.PD + sum(config.PD*config.g./(ones(config.Np_total,1)*(config.lambda+config.C_k)),2))'.*config.w;
            
            N_k = sum(config.w,2);
            round(N_k);

            [config.particles, config.w] = obj.resample(config.particles, (config.w/N_k)', config.resampling_strategy, config.Np);
            config.w = config.w'*N_k;
            obj.config = config;
        end
        
       % Resampling function
        function [xk, wk, idx] = resample(obj, xk, wk, resampling_strategy, Np_new)

            Np = length(wk);  % Np = number of particles

            % wk = wk./sum(wk); % normalize weight vector (already done)

            switch resampling_strategy
               case 'multinomial_resampling'
                  with_replacement = true;
                  idx = randsample(1:Np, Np, with_replacement, wk);
                %{
                  THIS IS EQUIVALENT TO:
                  edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
                  edges(end) = 1;                 % get the upper edge exact
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, idx] = histc(sort(rand(Np,1)), edges);
                %}
               case 'systematic_resampling'
                  % this is performing latin hypercube sampling on wk
                  edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
                  edges(end) = 1;                 % get the upper edge exact
                  u1 = rand/Np_new;
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, ~, idx] = histcounts(u1:1/Np_new:1, edges);
               otherwise
                  error('Resampling strategy not implemented\n')
            end
            xk = xk(:,idx);                    % extract new particles
            wk = repmat(1/Np_new, 1, Np_new)';          % now all particles have the same weight
        end
    end
    
    methods (Static)
        function D=mahalDist(x, m, C, use_log)
        % p=gaussian_prob(x, m, C, use_log)
        %
        % Evaluate the multi-variate density with mean vector m and covariance
        % matrix C for the input vector x.
        % Vectorized version: Here X is a matrix of column vectors, and p is 
        % a vector of probabilities for each vector.

            if nargin<4, use_log = 0; end

            d   = length(m);

            if size(x,1)~=d
               x=x';
            end
            N       = size(x,2);

            m       = m(:);
            M       = m*ones(1,N);
            denom   = (2*pi)^(d/2)*sqrt(abs(det(C)));
            invC    = inv(C);
            mahal   = sum(((x-M)'*invC).*(x-M)',2);   % Chris Bregler's trick

            switch use_log,
            case 2,
              D     = mahal;
            case 1,
              D     = -0.5*mahal - log(denom);
            case 0,
              numer = exp(-0.5*mahal);
              D     = numer/denom;
            otherwise
                error('Unsupported log type')
            end
        end
    end
end