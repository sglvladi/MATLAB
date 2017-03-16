classdef ParticleFilterMin2
% =====================================================================================
% Parameters:
% pf: structure with the following fields
%       
%       * Variables
%       -------------------
%       .Np               = number of particles
%       .w                = weights   (Np x T)
%       .particles        = particles (dim x Np x T)
%       .xhk              = estimated state
%       .z                = observation vector at time k (column vector)
%       .k                = iteration number
%       
%       * Handles
%       -------------------
%       .sys              = function handle to process equation
%       .obs              = function handle of the observation likelihood PDF p(y[k] | x[k])
%       .obs_model        = function handle to observation equation
%       .sys_noise        = function handle of a procedure that generates system noise
%       .gen_x0           = function handle of a procedure that samples from the initial pdf p_x0
%       .resample_strategy = resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
%

    properties
        pf
    end
    
    methods
        function obj = ParticleFilterMin2(prop)
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
                    %fprintf('Their std is %d \n', std(prop.particles(:,:),1));
%                     prop.particles(:,:)
%                     prop.w
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
             
            
            % Validate .obs
            if ~isfield(prop,'obs')
                error('Function handle for observation likelihood PDF p(y[k]|x[k]) (obs) is not given! Aborting...\n');
            end
            
            % Validate .z
            if ~isfield(prop,'z')
                fprintf('No initial observation supplied... Assuming "z = 0"\n');
                prop.z = 0;
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
            
            % Validate .clutter_flag
            if ~isfield(prop,'clutter_flag')
                fprintf('Function handle to process equation missing... Assuming "sys = .(x)x"..\n');
                prop.clutter_flag = 0;
            end
            
            % Validate .multi_flag
            if ~isfield(prop,'multi_flag')
                fprintf('Function handle to process equation missing... Assuming "sys = .(x)x"..\n');
                prop.multi_flag = 0;
            end
            
            if isfield(prop,'kf')
                fprintf('KF provided... Switching to Optimal proposal mode.');
                prop.optimal_prop_flag=1;
            else
                prop.optimal_prop_flag=0;
            end
            
            obj.pf = prop;
      
        end
        
        function pf = Predict(obj, pf)
             nx = size(pf.particles,1);               % number of states
             ny = size(pf.z,1);
             k = pf.k;
%             if k == 1
%                error('error: k must be an integer greater or equal than 2');
%             end
            
            if(pf.optimal_prop_flag)
                 pf.kf.s.z = pf.z;
                 pf.kf.s.x = (mean(pf.particles,2));
                 pf.kf.s.P = weightedcov(pf.particles',pf.w');
                 pf.kf.s = pf.kf.Iterate(pf.kf.s);
                 for i = 1:pf.Np
                    pf.particles(:,i) = mvnrnd(pf.kf.s.x,pf.kf.s.P); % Sample from optimal proposal
                 end
            else
                 %for i = 1:pf.Np
                 pf.particles(:,:) = pf.sys(k, pf.particles(:,:), pf.sys_noise(pf.particles(:,:))); % Simply propagate all particles
                 %end
            end
            
            if(pf.clutter_flag)
                % Expected likelihood variables 
                GateLevel   = 9;
                PG          = 0.989;      % probability of Gating
                PD          = 0.8;      % probability of Detection
                PointNum = size(pf.z,2); % number of measurements
                ObsDim = size(pf.z,1); % measurement dimensions
                C   = pi; % volume of the 2-dimensional unit hypersphere    
                
                trans_parts = pf.obs_model(pf.particles);   % Transform particles to observation space
                trans_parts = trans_parts + mvnrnd(zeros(size(trans_parts')),pf.R)';
                trans_mean = mean(trans_parts,2);
                
                pf.S = [std(trans_parts(1,:),pf.w')^2,0;0,std(trans_parts(2,:),pf.w')^2];
                pf.V_k = C*GateLevel^(ObsDim/2)*det(pf.S)^(1/2);   % volume of the validation region 
                
                % Particle-by-Particle validation matrix computation
%                Validation_matrix = zeros(1,PointNum);
%                 for i=1:pf.Np
%                     trans_part = trans_parts(:,i);
%                     
%                     for j=1:PointNum
%                         % distance
%                         if(Validation_matrix(1,j)==0)
%                             DistM(1,j)  = mahalDist(pf.z(:,j), mean(trans_part,2), S, 2);
%                             Validation_matrix(1,j)= DistM(1,j) < GateLevel;
%                         end
%                     end
%                     DistM2  = sum(((pf.z-mean(trans_part,2)*ones(1,PointNum)).*(pf.z-mean(trans_part,2)*ones(1,PointNum)))'/S,2)';
%                     Validation_matrix = Validation_matrix + DistM2 < GateLevel;
%                 end
                DistM  = sum(((pf.z-trans_mean*ones(1,PointNum)).*(pf.z-trans_mean*ones(1,PointNum)))'/pf.S,2)';

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % thresholding/ gating
                pf.Validation_matrix = DistM < GateLevel;
                if(isempty(find(pf.Validation_matrix,1)))
                    a=5;
                end
                pf.z_pred = mean(trans_parts,2);
            end
        end
        
        function pf = Update(obj, pf)
            if(pf.clutter_flag)
                % Expected likelihood variables 
                GateLevel   = 10;
                PG          = 0.989;      % probability of Gating
                PD          = 0.8;      % probability of Detection
                PointNum    = size(pf.z,2); % number of measurements
                ObsDim      = size(pf.z,1); % measurement dimensions
                C           = pi; % volume of the 2-dimensional unit hypersphere
                betta       = zeros(1,PointNum+1);  
                
                % Get valid measurement data indices
                DataInd = find(pf.Validation_matrix(1,:));
                
                % extract measurements
                z = pf.z(:,DataInd);
                [ObsDim,ValidDataPointNum] = size(z);
                % V_k = 10^2*pi;   % volume of the validation region: 2-D area of circle with radius=5
                lambda = ValidDataPointNum/pf.V_k; % expected rate of clutter per unit region (2 over entire volume)
                C1 = PD*PG*poisspdf(ValidDataPointNum-1,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum-1);
                C2 = (1-PD*PG)*poisspdf(ValidDataPointNum,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum);
                C11 = C1/(C1+C2);
                C22 = C2/(C1+C2);
            end
            % Initialize variables
            Np = pf.Np;                              % number of particles
            nx = size(pf.particles,1);               % number of states
            
            wk   = zeros(size(pf.w(:)));     % = zeros(Np,1);
            clutter_flag = pf.clutter_flag;
            % Update weights
            
            %for i = 1:Np                        
               
            if(~clutter_flag)    % If no clutter is present
                % Calculate new weights according to observation likelihood
                wk = pf.w .* pf.obs(k, pf.z, pf.particles(:,:));
            else
                z_pred = pf.obs_model(pf.particles(:,:));
                %Li_i = mvnpdf(z(:,:)', z_pred', pf.R);
                Li = zeros(size(z_pred,2),size(z, 2)+1);
                Li(:,1) = C22/(pf.V_k^(ValidDataPointNum));
                for i = 1:size(z, 2);
                    Li(:,i+1) = mvnpdf(z_pred', z(:,i)', pf.R)*C11/(pf.V_k^(ValidDataPointNum-1)*ValidDataPointNum);
                end

                % Calculate new weights according to expected likelihood
                wk = pf.w .* sum(Li(:,:),2); 

            end  
            %end;
            % Normalize weight vector
            pf.w = wk./sum(wk);
            
            % Calculate effective sample size: eq 48, Ref 1
            Neff = 1/sum(pf.w.^2);
            
            % Resampling
%             resample_percentage = 0.50;
%             Nt = resample_percentage*Np;
%             if Neff < Nt
               disp('Resampling ...')
               [pf.particles, pf.w] = obj.resample(pf.particles, pf.w, pf.resampling_strategy);
% %               {xk, pf.w} is an approximate discrete representation of p(x_k | y_{1:k})
%             end
            
            % Compute estimated state
            pf.xhk = zeros(nx,1);
            %for i = 1:Np
            pf.xhk(:,1) = sum(bsxfun(@times, pf.w(1,:)', pf.particles(:,:)),2);             
            
        end
        
        function pf = PredictMulti(obj, pf)
             Pdeath =  0.1;
             pf.ExistProb = (1 - Pdeath)*pf.ExistProb;
             nx = size(pf.particles,1);               % number of states
             ny = size(pf.z,1);
             k = pf.k;
             
%             if k == 1
%                error('error: k must be an integer greater or equal than 2');
%             end
            
            if(pf.optimal_prop_flag)
                 pf.kf.s.z = pf.z;
                 pf.kf.s.x = (mean(pf.particles,2));
                 pf.kf.s.P = weightedcov(pf.particles',pf.w');
                 pf.kf.s = pf.kf.Iterate(pf.kf.s);
                 for i = 1:pf.Np
                    pf.particles(:,i) = mvnrnd(pf.kf.s.x,pf.kf.s.P); % Sample from optimal proposal
                 end
            else
                 %for i = 1:pf.Np
                 pf.particles(:,:) = pf.sys(k, pf.particles(:,:), pf.sys_noise(pf.particles(:,:))); % Simply propagate all particles
                 %end
            end
            
            if(pf.clutter_flag)
                % Expected likelihood variables 
                GateLevel   = 9;
                PG          = 0.989;      % probability of Gating
                PD          = 0.8;      % probability of Detection
                PointNum = size(pf.z,2); % number of measurements
                ObsDim = size(pf.z,1); % measurement dimensions
                C   = pi; % volume of the 2-dimensional unit hypersphere    
                
                trans_parts = pf.obs_model(pf.particles);   % Transform particles to observation space
                trans_parts = trans_parts + mvnrnd(zeros(size(trans_parts')),pf.R)';
                trans_mean = mean(trans_parts,2);
                
                %pf.S = [std(trans_parts(1,:),pf.w')^2,0;0,std(trans_parts(2,:),pf.w')^2];
                pf.S = weightedcov(trans_parts',pf.w');
                pf.V_k = C*GateLevel^(ObsDim/2)*det(pf.S)^(1/2);   % volume of the validation region 
                
                % Particle-by-Particle validation matrix computation
%                Validation_matrix = zeros(1,PointNum);
%                 for i=1:pf.Np
%                     trans_part = trans_parts(:,i);
%                     
%                     for j=1:PointNum
%                         % distance
%                         if(Validation_matrix(1,j)==0)
%                             DistM(1,j)  = mahalDist(pf.z(:,j), mean(trans_part,2), S, 2);
%                             Validation_matrix(1,j)= DistM(1,j) < GateLevel;
%                         end
%                     end
%                     DistM2  = sum(((pf.z-mean(trans_part,2)*ones(1,PointNum)).*(pf.z-mean(trans_part,2)*ones(1,PointNum)))'/S,2)';
%                     Validation_matrix = Validation_matrix + DistM2 < GateLevel;
%                 end
                DistM  = sum(((pf.z-trans_mean*ones(1,PointNum)).*(pf.z-trans_mean*ones(1,PointNum)))'/pf.S,2)';

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % thresholding/ gating
                pf.Validation_matrix = DistM < GateLevel;
                if(isempty(find(pf.Validation_matrix,1)))
                    a=5;
                end
                % Get valid measurement data indices
                DataInd = find(pf.Validation_matrix(1,:));

                % extract measurements
                z = pf.z(:,DataInd);
                [ObsDim,ValidDataPointNum] = size(z);
                
                lambda = ValidDataPointNum/pf.V_k; % expected rate of clutter per unit region (2 over entire volume)
                C1 = PD*PG*poisspdf(ValidDataPointNum-1,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum-1);
                C2 = (1-PD*PG)*poisspdf(ValidDataPointNum,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum);
                C11 = C1/(C1+C2);
                C22 = C2/(C1+C2);
                
                pf.z_pred = mean(trans_parts,2);
                %% Compute Association Likelihoods 
                z_pred = pf.obs_model(pf.particles(:,:));
                pf.Li = zeros(pf.Np, ValidDataPointNum);
                try
                    for i = 1:ValidDataPointNum;
                        pf.Li(:,i) = mvnpdf(z_pred', pf.z(:,i)', pf.R);
                    end
                    pf.Li = sum(pf.Li, 1)/pf.Np;
                catch
                    disp('Association Likelihood error!!');
                end
            end
        end
        
        function pf = UpdateMulti(obj, pf)
            PG          = 0.989;      % probability of Gating
            PD          = 0.8;      % probability of Detection
            % Get valid measurement data indices
            DataInd = find(pf.Validation_matrix(1,:));

            % extract measurements
            z = pf.z(:,DataInd);
            [ObsDim,ValidDataPointNum] = size(z);
            
            lambda = ValidDataPointNum/pf.V_k; % expected rate of clutter per unit region (2 over entire volume)
            C1 = PD*PG*poisspdf(ValidDataPointNum-1,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum-1);
            C2 = (1-PD*PG)*poisspdf(ValidDataPointNum,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum);
            C11 = C1/(C1+C2);
            C22 = C2/(C1+C2);
          
            % Initialize variables
            Np = pf.Np;                              % number of particles
            nx = size(pf.particles,1);               % number of states
            
            wk   = zeros(size(pf.w(:)));     % = zeros(Np,1);
            clutter_flag = pf.clutter_flag;
            
            %% Compute Association Likelihoods 
            z_pred = pf.obs_model(pf.particles(:,:));
            Li = zeros(size(z_pred,2),ValidDataPointNum+1);
            try
                Li(:,1) = ones(size(z_pred,2),1)*pf.betta(1);%*C22/(pf.V_k^(ValidDataPointNum));
            catch
                disp('error');
            end
            try
                if(size(pf.betta,2)~=1)
                    for i = 1:size(z, 2);
                        Li(:,i+1) = mvnpdf(z_pred', z(:,i)', pf.R)*pf.betta(i+1);%*C11/(pf.V_k^(ValidDataPointNum-1)*ValidDataPointNum);
                    end
                end
            catch
                disp('Association Likelihood error!!');
            end
            
            % Calculate new weights according to expected likelihood
            wk = pf.w .* sum(Li(:,:),2); 
            
            % Normalize weight vector
            pf.w = wk./sum(wk);
            
            % Calculate effective sample size: eq 48, Ref 1
            Neff = 1/sum(pf.w.^2);
            
            % Resampling
%             resample_percentage = 0.50;
%             Nt = resample_percentage*Np;
%             if Neff < Nt
               disp('Resampling ...')
               [pf.particles, pf.w] = obj.resample(pf.particles, pf.w, pf.resampling_strategy);
% %               {xk, pf.w} is an approximate discrete representation of p(x_k | y_{1:k})
%             end
            
            % Compute estimated state
            pf.xhk(:,1) = sum(bsxfun(@times, pf.w(1,:)', pf.particles(:,:)),2);
            
            % Update existence probability
            if(pf.clutter_flag)
                pf.ExistProb = (1-pf.ExistProb)*sum(Li(:,2:end))/(sum(Li(:,1))+sum(Li(:,2:end))) + pf.ExistProb; 
                disp(pf.ExistProb);
            end
                
        end
        
        function pf = PredictSearch(obj, pf)
            Pbirth = 0.005;
            Pdeath =  0.1;
            nx = size(pf.particles,1);               % number of states
            ny = size(pf.z,1);
            k = pf.k;
%             if k == 1
%                error('error: k must be an integer greater or equal than 2');
%             end
            pf.ExistProb =  Pbirth*(1-pf.ExistProb) + (1-Pdeath)*pf.ExistProb;
            
            if(pf.optimal_prop_flag)
                 pf.kf.s.z = pf.z;
                 pf.kf.s.x = (mean(pf.particles,2));
                 pf.kf.s.P = weightedcov(pf.particles',pf.w');
                 pf.kf.s = pf.kf.Iterate(pf.kf.s);
                 for i = 1:pf.Np
                    pf.particles(:,i) = mvnrnd(pf.kf.s.x,pf.kf.s.P); % Sample from optimal proposal
                 end
            else
                 %for i = 1:pf.Np
                 pf.particles(:,:) = pf.sys(k, pf.particles(:,:), pf.sys_noise(pf.particles(:,:))); % Simply propagate all particles
                 %end
            end
            
            if(pf.clutter_flag)
                % Expected likelihood variables 
                GateLevel   = 9;
                PG          = 0.989;      % probability of Gating
                PD          = 0.8;      % probability of Detection
                PointNum = size(pf.z,2); % number of measurements
                ObsDim = size(pf.z,1); % measurement dimensions
                C   = pi; % volume of the 2-dimensional unit hypersphere    
                
                trans_parts = pf.obs_model(pf.particles);   % Transform particles to observation space
                trans_parts = trans_parts + mvnrnd(zeros(size(trans_parts')),pf.R)';
                trans_mean = mean(trans_parts,2);
                
                pf.S = [std(trans_parts(1,:),pf.w')^2,0;0,std(trans_parts(2,:),pf.w')^2];
                pf.V_k = C*GateLevel^(ObsDim/2)*det(pf.S)^(1/2);   % volume of the validation region 
                
                % Particle-by-Particle validation matrix computation
%                Validation_matrix = zeros(1,PointNum);
%                 for i=1:pf.Np
%                     trans_part = trans_parts(:,i);
%                     
%                     for j=1:PointNum
%                         % distance
%                         if(Validation_matrix(1,j)==0)
%                             DistM(1,j)  = mahalDist(pf.z(:,j), mean(trans_part,2), S, 2);
%                             Validation_matrix(1,j)= DistM(1,j) < GateLevel;
%                         end
%                     end
%                     DistM2  = sum(((pf.z-mean(trans_part,2)*ones(1,PointNum)).*(pf.z-mean(trans_part,2)*ones(1,PointNum)))'/S,2)';
%                     Validation_matrix = Validation_matrix + DistM2 < GateLevel;
%                 end
                DistM  = sum(((pf.z-trans_mean*ones(1,PointNum)).*(pf.z-trans_mean*ones(1,PointNum)))'/pf.S,2)';

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % thresholding/ gating
                pf.Validation_matrix = DistM < GateLevel;
                if(isempty(find(pf.Validation_matrix,1)))
                    a=5;
                end
                pf.z_pred = mean(trans_parts,2);
            end
        end
        
        function pf = UpdateSearch(obj, pf)
            if(pf.clutter_flag)
                % Expected likelihood variables 
                GateLevel   = 10;
                PG          = 0.989;      % probability of Gating
                PD          = 0.8;      % probability of Detection
                PointNum    = size(pf.z,2); % number of measurements
                ObsDim      = size(pf.z,1); % measurement dimensions
                C           = pi; % volume of the 2-dimensional unit hypersphere
                betta       = zeros(1,PointNum+1);  
                
                % Get valid measurement data indices
                DataInd = find(pf.Validation_matrix(1,:));
                
                % extract measurements
                z = pf.z(:,DataInd);
                [ObsDim,ValidDataPointNum] = size(z);
                % V_k = 10^2*pi;   % volume of the validation region: 2-D area of circle with radius=5
                lambda = ValidDataPointNum/pf.V_k; % expected rate of clutter per unit region (2 over entire volume)
                C1 = PD*PG*poisspdf(ValidDataPointNum-1,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum-1);
                C2 = (1-PD*PG)*poisspdf(ValidDataPointNum,lambda*pf.V_k);%/pf.V_k^(ValidDataPointNum);
                C11 = C1/(C1+C2);
                C22 = C2/(C1+C2);
            end
            % Initialize variables
            Np = pf.Np;                              % number of particles
            nx = size(pf.particles,1);               % number of states
            
            wk   = zeros(size(pf.w(:)));     % = zeros(Np,1);
            clutter_flag = pf.clutter_flag;
            
            % Update weights       
            if(~clutter_flag)    % If no clutter is present
                % Calculate new weights according to observation likelihood
                wk = pf.w .* pf.obs(k, pf.z, pf.particles(:,:));
            else
                z_pred = pf.obs_model(pf.particles(:,:));
                %Li_i = mvnpdf(z(:,:)', z_pred', pf.R);
                Li = zeros(size(z_pred,2),size(z, 2)+1);
                Li(:,1) = C22/(pf.V_k^(ValidDataPointNum));
                for i = 1:size(z, 2);
                    Li(:,i+1) = mvnpdf(z_pred', z(:,i)', pf.R)*C11/(pf.V_k^(ValidDataPointNum-1)*ValidDataPointNum);
                end

                % Calculate new weights according to expected likelihood
                if(pf.ExistProb<0.9)
                    wk = pf.w .* sum(Li(:,:),2);
                    pf.ExistProb = (1-pf.ExistProb)*sum(Li(:,2:end))/(sum(Li(:,1))+sum(Li(:,2:end))) + pf.ExistProb; 
                    disp(pf.ExistProb);
                else
                    % When promoting a track, condition on the most likely
                    % measurement to avoid multi-modality
                    [betta_max, betta_max_ind]=max(betta);
                    wk = pf.w .* Li(:,betta_max_ind);
                    pf.ExistProb = sum(Li(:,2:end))/(sum(Li(:,1))+sum(Li(:,2:end))); 
                    disp(pf.ExistProb);
                end

            end  
            %end;
            % Normalize weight vector
            pf.w = wk./sum(wk);
            
            % Calculate effective sample size: eq 48, Ref 1
            Neff = 1/sum(pf.w.^2);
            
            % Resampling
%             resample_percentage = 0.50;
%             Nt = resample_percentage*Np;
%             if Neff < Nt
               disp('Resampling ...')
               [pf.particles, pf.w] = obj.resample(pf.particles, pf.w, pf.resampling_strategy);
% %               {xk, pf.w} is an approximate discrete representation of p(x_k | y_{1:k})
%             end
            
            % Compute estimated state
            pf.xhk = zeros(nx,1);
            %for i = 1:Np
            pf.xhk(:,1) = sum(bsxfun(@times, pf.w(1,:)', pf.particles(:,:)),2);             
            
        end
        
        function pf = Iterate(obj, pf)
            k = pf.k;
            if k == 1
               error('error: k must be an integer greater or equal than 2');
            end
            
            if(pf.optimal_prop_flag)
                 pf.kf.s.z = pf.z;
                 pf.kf.s.x = (mean(pf.particles,2));
                 pf.kf.s.P = weightedcov(pf.particles',pf.w');
                 pf.kf.s = pf.kf.Iterate(pf.kf.s);
                 %for i = 1:pf.Np
                    pf.particles(:,:) = mvnrnd((pf.kf.s.x*ones(1, pf.Np))',pf.kf.s.P)'; % Sample from optimal proposal
                 %end
            else
                 %for i = 1:pf.Np
                 pf.particles(:,:) = pf.sys(k, pf.particles(:,:), pf.sys_noise(pf.particles(:,:))); % Simply propagate all particles
                 %end
            end
            
            % Initialize variables
            Np = pf.Np;                              % number of particles
            nx = size(pf.particles,1);               % number of states
            
            if(~pf.clutter_flag)    % If no clutter is present
                % Calculate new weights according to observation likelihood
                pf.w = pf.w .* pf.obs(k, pf.z, pf.particles(:,:));
             else
%                  Li = zeros(ValidDataPointNum+1, 1);
%                  z_pred = pf.obs_model(pf.particles(:,i));
% 
%                 %% Compute Association Likelihoods 
%                 Li(ValidDataPointNum+1,1) = C22;
%                 Li(1:ValidDataPointNum) = C11*mvnpdf(z(:,:)', z_pred', pf.R);
%                 % Calculate new weights according to expected likelihood
%                 pf.w(i) = pf.w(i) * sum(Li,1); %obj.expectedLikelihood(pf.z, z_pred, pf.R);
            end  
            % Normalize weight vector
            pf.w = pf.w./sum(pf.w);
            
            % Calculate effective sample size: eq 48, Ref 1
            Neff = 1/sum(pf.w.^2);
            
            % Resampling
             resample_percentage = 0.50;
             Nt = resample_percentage*Np;
            if Neff < Nt
               disp('Resampling ...')
               [pf.particles, pf.w] = obj.resample(pf.particles, pf.w, pf.resampling_strategy);
% %               {xk, pf.w} is an approximate discrete representation of p(x_k | y_{1:k})
             end
            
            % Compute estimated state
            pf.xhk = zeros(nx,1);
            % Compute estimated state
            pf.xhk(:,1) = sum(bsxfun(@times, pf.w(:,1)', pf.particles(:,:)),2);
            
        end
        
        
        % Resampling function
        function [xk, wk, idx] = resample(obj, xk, wk, resampling_strategy)

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
                  u1 = rand/Np;
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, ~, idx] = histcounts(u1:1/Np:1, edges);
               otherwise
                  error('Resampling strategy not implemented\n')
            end;
            xk = xk(:,idx);                    % extract new particles
            wk = repmat(1/Np, 1, Np)';          % now all particles have the same weight
        end
        
        function expLi = expectedLikelihood(obj,DataList, z_pred, R)
            GateLevel   = 5;
            PG          = 0.918;      % probability of Gating
            PD          = 0.8;      % probability of Detection
            PointNum = size(DataList,2); % number of measurements
            ObsDim = size(DataList,1); % measurement dimensions 
            V_k = 5^2*pi;   % volume of the validation region: 2-D area of circle with radius=5
            lambda = 1/V_k; % expected rate of clutter per unit region (2 over entire volume)

            %% Compute Association Likelihoods 
            Li = zeros(PointNum+1, 1);
            Li(PointNum+1,1) = (1-PD*PG)*poisspdf(PointNum,lambda*V_k)/V_k^(PointNum);
            %Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG);
            Li(1:PointNum) = (1/PointNum)*PD*PG*poisspdf(PointNum-1,lambda*V_k)*mvnpdf(DataList(:,:)', z_pred', R)/(PG*V_k^(PointNum-1));
%             parfor i=1:PointNum
%                 z = DataList(:,i);
%                 Li(i,1) = (1/PointNum)*PD*PG*poisspdf(PointNum-1,lambda*V_k)*mvnpdf(z, z_pred, R)/(PG*V_k^(PointNum-1));
%             end
            
            expLi = sum(Li,1);
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