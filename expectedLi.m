function expLi = expectedLi(DataList, z_pred, S, R)
    GateLevel   = 5;
    PG          = 0.918;      % probability of Gating
    PD          = 0.8;      % probability of Detection
    PointNum = size(DataList,2); % number of measurements
    ObsDim = size(DataList,1); % measurement dimensions
    C   = pi; % volume of the 2-dimensional unit hypersphere     
    V_k = C*GateLevel^(ObsDim/2)*det(S)^(1/2);   % volume of the validation region 
    
    %% Possibly perform gating here!
    
    %% Compute Association Likelihoods 
    Li = zeros(PointNum, 1);
    %Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG);
    for i=1:PointNum
        z = DataList(:,i);
        Li(i,1) = mvnpdf(z, z_pred, S)*PD/(PointNum/V_k);
    end
    
    %% Compute Observation Likelihoods 
    Li_k = zeros(PointNum, 1);
    %Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG);
    for i=1:PointNum
        z = DataList(:,i);
        Li_k(i,1) = mvnpdf(z, z_pred, R);
    end
    
    % Compute association probabilities
    betta(1:PointNum) = Li(1:PointNum)./(1-PG*PD+sum(Li,1));
    betta(PointNum+1) = (1-PG*PD)/(1-PG*PD+sum(Li,1));
    expLi = betta(PointNum+1)/(V_k^PointNum) + sum(betta(1:PointNum)'.*Li_k(:,1),1)/(PG*PD*V_k^(PointNum-1));
end