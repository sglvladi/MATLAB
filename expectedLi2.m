function expLi = expectedLi2(DataList, z_pred, S, R)
    GateLevel   = 5;
    PG          = 0.918;      % probability of Gating
    PD          = 0.8;      % probability of Detection
    PointNum = size(DataList,2); % number of measurements
    ObsDim = size(DataList,1); % measurement dimensions
    C   = pi; % volume of the 2-dimensional unit hypersphere     
    V_k = C*GateLevel^(ObsDim/2)*det(S)^(1/2);   % volume of the validation region 
    
    %% Compute Association Likelihoods 
    Li = zeros(PointNum+1, 1);
    %Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG);
    for i=1:PointNum
        z = DataList(:,i);
        Li(i,1) = PointNum*mvnpdf(z, z_pred, R)*PD*PG/(V_k^(PointNum-1));
    end
    Li(PointNum+1,1) = (1-PD*PG)/(V_k^(PointNum));
    
    expLi = sum(Li.^2)/sum(Li)
end