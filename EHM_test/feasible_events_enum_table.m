% 2 targets, 1 common measurement
% Omega = [...
%     1 1 1 ];

% 2 targets, 2 common measurement
Omega = [...
    1 1 
    1 1];

% Example from "Fast Mutual Exlusion", S. Maskell et. al., 2004
% Omega = [...
%     1 1 0 0 0
%     1 1 1 0 0
%     1 1 1 1 0
%     1 0 0 1 1];

% Omega = [...
%     1 0 1 0  
%     1 1 1 0 
%     1 1 0 0 
%     1 1 0 1 
%     1 0 0 1 ];
%  

%  Omega = [...
%     1 1 0 1  
%     1 0 1 0 
%     1 1 0 0 
%     1 0 0 1 
%     1 1 1 0 ];
% 
% Omega = [...
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1];


% Omega = [...
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1];

% Omega = [...
%     1 1 1 1;
%     1 0 1 1;
%     1 0 0 1;
%     1 0 0 0];

% Example 2.1
% Omega = [...
%     1 1 0 0 1 0;
%     1 1 1 1 0 1;
%     1 1 1 0 0 0;
%     1 0 0 1 1 1];

%Example 6
% Omega = [...
%     1 1 1 1 0 0 0 0 1 0 1;
%     1 0 1 1 0 0 0 1 1 0 0;
%     1 0 0 0 1 1 1 0 0 0 0;
%     1 0 0 1 1 1 1 0 0 1 1;
%     1 0 0 1 1 0 0 0 0 1 1];

% Example 7
% Omega = [...
%     1 1 1 1 0 0 0 0 1 0 1;
%     1 0 1 1 0 0 0 1 1 0 0;
%     1 0 0 0 1 1 1 0 0 0 0;
%     1 0 0 1 1 1 1 0 0 1 1;
%     1 0 0 1 1 0 0 0 0 1 1;
%     1 0 0 0 1 0 0 1 0 0 0;
%     1 0 0 0 1 1 1 1 1 0 0];

% Example 8
% Omega = [...
%     1 1 1 1 0 1 0 0 1 0;
%     1 0 1 1 0 0 0 0 0 1;
%     1 0 0 0 1 1 0 0 0 0;
%     1 1 0 0 0 0 1 0 1 1;
%     1 0 0 0 0 0 0 1 0 0;
%     1 0 0 0 0 0 0 1 1 1;
%     1 0 1 1 1 1 0 0 0 0];

% Example 8.1
% Omega = [...
%     1 1 1 1 1 1 0 1 0 0;
%     1 0 1 1 0 0 0 0 0 0;
%     1 0 0 0 1 0 1 0 0 0;
%     1 0 0 0 0 0 1 0 1 1;
%     1 0 0 0 0 0 0 1 0 0;
%     1 0 0 0 0 0 0 1 1 1;
%     1 0 1 1 0 0 0 0 0 0];

% Example 8.2
% Omega = [...
%     1 1 1 1 1 1 0 0 0 0;
%     1 0 1 1 0 0 0 0 0 0;
%     1 0 0 0 0 0 1 1 0 0;
%     1 0 0 0 0 0 0 0 1 1;
%     1 0 0 0 0 0 0 1 0 0;
%     1 0 0 0 0 0 0 0 1 0;
%     1 0 1 1 0 1 0 0 0 0];

% % Example 8.3
% Omega = [...
%     1 1 1 1 1 1 0 0 0 0;
%     1 0 1 1 0 0 0 0 0 0;
%     1 0 0 0 0 0 1 1 0 1;
%     1 0 0 0 0 0 0 0 1 1;
%     1 0 0 0 0 0 0 1 0 0;
%     1 0 0 0 0 0 0 0 1 0;
%     1 0 1 1 0 1 0 0 0 0];

% Example 9
% Omega = [...
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1;
%     1 1 1 1 1 1 1 1 1 1];
[m,n] = size(Omega);
Nt = n-1;

%% Number of validated measurements (sum over targets - columns)
indvm = (sum(Omega,2) > 0)';
mk = sum(indvm);
Omegaf = zeros(mk,Nt+1);
Omegaf(1:mk,:) = Omega(indvm,:);

%% Generation of validation matrix is in Ref. 1, p. 390

% Combinations of unique (feasible) associations for all tracks (all rows)
sc = sum(Omegaf,2);
% All possible combinations without the constraint of at most one
% measurement originaged from a target
c = prod(sc);

% Association events
theta = cell(c,1);

% Generation of feasible events
ind = cell(mk,1);
Om = zeros(c*mk,size(Omegaf,2));
Oms = zeros(c*mk,size(Omegaf,2));

delta = zeros(1,c*Nt);
tau = zeros(c*mk,1);
phi = zeros(c,1);

% Indices of tracks
for j = 1:mk
    ind{j} = find(Omegaf(j,:) == 1);
end

fprintf('\nRunning MO JPDA with enumerating table...\n===============================================>');
tic;
%% JPDA with enumerating table
% Matrix with all association events that satisfy:
% one source for each measurement (one target market in each row)
for j = 1:mk % 1 to 4
    nt = sc(j); % 4 | 3 | 3 | 2
    nl = prod(sc(1:j)); % 4 | 4*3 (=12) | 4*3*3 (=36) | 4*3*3*2 (=72)
    pf = c*mk/nl;     % 72c*4m/4t (=72) | 72c*4m/(4*3)t (=24) 
    % 72c*4m/(4*3*3)t (=8) | 72c*4m/(4*3*3*2)t (=4)
    base = repmat(eye(nt),nl/nt,1); % eye(4) | eye(3) x 4 | eye(3) x 4*3 | eye(2) x 4*3*3
    for k = 1:nl % 1 to 4 | 1 to 4*3 | 1 to 4*3*3 | 1 to 4*3*3*2
        for L = 1:pf/mk
            Om(pf*k -pf +j +mk*L -mk,ind{j}) = base(k,1:nt);
        end
    end
end

% Eliminate combinations with more than one
% measurement originaged from a target
k = 0;
for i = 1:c
    Omi = Om(i*mk -mk +1:i*mk,:);
    deltai = sum(Omi(:,2:end),1); % Defined for t = 1..Nt (exclude target 0)
    if ~(deltai > 1)
        k = k+1;
        Oms(k*mk -mk +1:k*mk,:) = Omi;
        tau(k*mk -mk +1:k*mk,1) = sum(Omi(:,2:end),2); % Defined for t = 1..Nt (exclude target 0)
        delta(1,k*Nt -Nt +1:k*Nt) = deltai;
        phi(k,1) = sum((ones(mk,1)-tau(k*mk -mk +1:k*mk,1)),1);
        
        theta{k,1} = Omi;
    end
end

% Allocate memory
Omf = zeros(k*mk,size(Omegaf,2));
deltaf = zeros(1,k*Nt);
tauf = zeros(k*mk,1);
phif = zeros(k,1);

thetaf = cell(k,1);

% Transfer the feasible events
Omf(1:k*mk,:) = Oms(1:k*mk,:);
deltaf(1,1:k*Nt) = delta(1,1:k*Nt);
tauf(1:k*mk,1) = tau(1:k*mk,1);
phif(1:k,1) = phi(1:k,1);
thetaf(1:k,1) = theta(1:k,1);
%% JPDA with table -- end
%% Define parameters
% Dimension of output vector
nz = 2;

% Volume of validation region
PDt = 0.99;
lambda = 0.3317;
gamma_ = chi2inv(0.99,nz);

%% Calculate marginal association probabilities for targets/tracks
%F = zeros(size(Omegaf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is used element-wise in the loop - by Flávio (13/07/2016)
% (no need to bother about the term lambda.(1-Pd) in this matrix)
%F(:,1) = zeros(mk,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = [0; 0];
S = eye(nz)/10;

cnz = pi^(nz/2)/gamma(nz/2 + 1);
V(1,1) = 0;
V(1,2:end) = cnz*sqrt(det(gamma_*S))*ones(size(V(1,2:end)));

% This is fictitious (from line 194 to 204) - for testing purposes
% for t = 2:Nt+1
% 	for j = 1:mk
% 		if Omegaf(j,t) == 1
% 			% F(j,t) = mvcpdf(zkv{j}, zhkp{t}, cov_zhkp{t});
%             z = (1/(j*t))*ones(nz,1);
%             F(j,t) = mvnpdf(z', mu', S);
%             %F(j,t) = Omegaf(j,1);
% 		end
% 	end
% end

beta = zeros(size(Omegaf));

%Target t = 0 (no detection) shall not be included
%Already taken into account considering the probability of detection
for t = 1:Nt+1
    for j = 1:mk
        betatj = 0;
        % ctj = 0;
        for i = 1:k
            P1 = 1; P2 = 1;
            % Just perform the calculation if it is part of the event to
            % compute marginal
            if Omf(i*mk -mk +j,t) == 1
                % For targets where
                ind = find(sum(Omf(i*mk -mk +1:i*mk,:),1) == 1);
                for ti = ind
                    a1 = ((lambda^-1)*F(:,ti)).^(tauf(i*mk -mk +1:i*mk,1).*Omf(i*mk -mk +1:i*mk,ti));
                    P1 = P1*prod(a1);
                end
                
                a2 = (PDt.^deltaf(1,i*Nt -Nt +1:i*Nt)).*...
                    ((1 -PDt).^(1 -deltaf(1,i*Nt -Nt +1:i*Nt)));
                P2 = prod(a2);
                
                P = P1*P2;
                % Include if the assignment is such that wjt = 1
                % if Omf(i*mk -mk +j,t) == 1
                % if Omegaf(j,t) == 1
                betatj = betatj + P;
                % end
                % Normalization (sum over all valid assignments)
                % ctj = ctj + P;
            end
        end
        beta(j,t) = betatj;
    end
end

% Normalise
for j = 1:mk
    beta(j,:) = beta(j,:)/sum(beta(j,:),2);
end
beta_mo = beta
toc;

% Compose the likelihood matrix as expected by your function - by Flávio (13/07/2016)
L = zeros(size(F));
L(:,1) = lambda*(1-PDt)*ones(mk,1);
L(:,2:Nt+1) = PDt*F(:,2:Nt+1);

% Transpose Omega
Omega_trans = [ones(size(Omega,2)-1,1), Omega(:, 2:end)'];

% Transpose likelihood matrix
L_trans = zeros(size(Omega_trans,1), size(Omega_trans,2));
L_trans(:,1) = lambda*(1-PDt)*ones(Nt,1);
L_trans(:,2:end) = L(:,2:end)';

% Generate results for Solution 1 (TO Tree based)
% ===============================================>
fprintf('\nRunning Solution 1 (TO Tree based)...\n===============================================>');
tic;
TreeObj = buildTree(Omega_trans, L_trans);
beta_to = TreeObj.betta
toc;
drawTreeObj(TreeObj, 'TO', 'Solution 1 Tree');

% Generate results for Solution 2 (MO-EHM)
% ===============================================>
fprintf('\nRunning Solution 2 (MO-EHM)...\n===============================================>');
tic;
NetObj = buildEHMnet2_flavio(Omega, L); 
beta_mo = NetObj.betta
beta_to = NetObj.betta_trans
toc;
p = drawNetObj(NetObj, 'MO', 'Solution 2 Net');

% Generate results for Solution 3
% ===============================================>
fprintf('\nRunning Solution 3 (TO-EHM)...\n===============================================>');
tic;
NetObj2 = buildEHMnet_trans(Omega_trans, L_trans);
beta_to = NetObj2.betta
toc;
D = drawNetObj(NetObj2, 'TO', 'Solution 3 Net');
;


        
