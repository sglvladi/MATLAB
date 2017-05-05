function pD = compute_pD(model,X)

% compute pD at a vector of values X
% this function is problem dependent
% here it is a Gaussian which depends only on position 

mid= [0; 0]; 
% cov= diag([4000,4000].^2);  % the test used before.
cov= diag([6000,6000].^2);

M= size(X,2);
P= model.pos*X(1:4,:);  %coordinate associated with X
e_sq= sum( (diag(1./diag(sqrt(cov)))*(P-repmat(mid,[1 M]))).^2 );
%pD= exp(-e_sq/2)';
pD= model.Pd*exp(-e_sq/2)'; %use this to scale the peak of the function i.e. average value by P_D)
% pD= model.Pd*ones(1,size(X,2))';        %use this for a constant P_D for comparison with EK and UK approximations
