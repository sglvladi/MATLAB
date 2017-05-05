
% An unified script for propagating particles according to the state
% procee function; for different type of target models and noises

if size(px,1)==4
    
%     if size(Noise,1)==2
%         for j=1:N
%             xPre(:,j) = F(px(:,j)) + G*(Noise.*randn(2,1)); % a faster way is don't use loop but deal with each dim separately
%             y(:,j) = H*xPre(:,j);
%         end
%     elseif size(Noise,1)==4
%         Mu =zeros(4,1);
%         for j=1:N
%             xPre(:,j) = F(px(:,j)) + mvnrnd(Mu,Noise)';
%             y(:,j) = H*xPre(:,j);
%         end
%     end

%  % the following is for fast computing: assume F = @(a)[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1]*a;
    T=1;
    if size(Noise,1)==2
        Pnoise=G*(repmat(Noise,1,N).*randn(2,N));
%         Pnoise=repmat((G*Noise),1,N).*randn(4,N);
    elseif size(Noise,1)==4
        Pnoise=mvnrnd(zeros(4,1),Noise,N)'; 
    end
     xPre(1,:)= px(1,:)+ px(2,:)*T + Pnoise(1,:);
     xPre(2,:)= px(2,:)+ Pnoise(2,:) ;
     xPre(3,:)= px(3,:)+ px(4,:)*T + Pnoise(3,:);
     xPre(4,:)= px(4,:)+ Pnoise(4,:);

elseif size(px,1)==5 
    
%     for j=1:N
%         xPre(:,j) = F(px(:,j)) + [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
%         y(:,j) = H*xPre(1:4,j);
%     end

    % % the following is for fast computing, assume
    % F = @(a)[1 sin(a(5,:)*T)/a(5,:)      0    -((1-cos(a(5,:)*T))/a(5,:)) 0
%          0 cos(a(5,:)*T)             0    -sin(a(5,:)*T)              0
%          0 (1-cos(a(5,:)*T))/a(5,:)  1    sin(a(5,:)*T)/a(5,:)        0
%          0 sin(a(5,:)*T)             0    cos(a(5,:)*T)               0
%          0       0                   0    0                           1]*a;
     T=1;
     if size(Noise.vel,1)==2
        Pnoise =  [G*(repmat(Noise.vel,1,N).*randn(2,N)); randn(1,N) * Noise.turn];
        xPre(5,:)=  px(5,:)+ Pnoise(5,:);
        px5 = px(5,:);   
     elseif size(Noise,1)==4
%         Pnoise =  sqrtm(Noise)*randn(size(Noise,2),1);
        Pnoise =  repmat(Noise,1,N)*randn(4,N);
        px5 = Noise_turn_rate; % need to specify the turn rate
     else
         disp('unstored state dynamics')
     end
     xPre(1,:)= px(1,:)+ (sin(px5*T)./px5).*px(2,:)-(((1-cos(px5*T))./px5)).*px(4,:)+ Pnoise(1,:);
     xPre(2,:)= cos(px5*T).*px(2,:) -sin(px5*T).*px(4,:)+ Pnoise(2,:);
     xPre(3,:)= ((1-cos(px5*T))/px5).*px(2,:)+px(3,:)+(sin(px5*T)/px5).*px(4,:)+ Pnoise(3,:);
     xPre(4,:)=  sin(px5*T).*px(2,:)+cos(px5*T).*px(4,:)+ Pnoise(4,:);

end

% called by CompWeight_

% 
% %% input;
% % Y:        Observations
% % w:        The weight of particles
% % x:        The state of previous particles
% % F,G,Q:    The simulation model parameters, the state equation is:
% %                  x_k=F*x_(k-1) + G*(Q.*randn(2,1))
% % H:        The observation model, the observation equation is;
% %                  y=H*x
% % birth:    new target birthing model
% % Msample:  The number of particles allocated to per new targets, 
% %           in our case we use Msample=Np    
% % Wd        the measurement noise
% % Ps:       The survival probability of taregts
% %% output;
% % xPre:     Predicted weight of particles,(new particles added)
% % glik:     The likelihood matrix, appearing as glik(j,z) 
% %           where j index particle, and z index  measurements. 
% %           glik(j,z) is corresponding to the g(z|j) in my paper;
% % wPre:     Predicted weight of particles, i.e. w_(k|k)
% 


