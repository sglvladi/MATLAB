
function [xPre,glik,wPre]=CompWeight_distance_bearing(Y,pw,px,H,F,G,Noise,sensor,r_noise,Birth,Np,Wd,Ps)

% % Computation the weight, based on active sensor (range-bearing observation)
% % main codes prepared for myphd.m

% %% input;
% % Y:        Observations
% % pw:       The weight of particles
% % px:       The state of previous particles
% % F:        The dynamic/process model   
% %                  x_k=F(x_(k-1)) 
% % G,q_noise: The process noise
% %                  + [G*(Q.*randn(2,1)); turn_noise]
% % H:        The observation model, the observation equation is;
% %                  y=H*x(1:4,:)
% % sensor:   the sensor localization
% % r_noise,  the obsevation noise
% % Birth:    new target Birthing model
% % Np:       The number of particles allocated to per new targets, 
% %           in our case we use Msample=Np    
% % Wd        the measurement noise
% % Ps:       The survival probability of taregts
% % output;
% % xPre:     Predicted weight of particles,(new particles added)
% % glik:     The likelihood matrix, appearing as glik(j,z) 
% %           where j index particle, and z index  measurements. 
% %           glik(j,z) is corresponding to the g(z|j) in my paper;
% % wPre:     Predicted weight of particles, i.e. w_(k|k)
%%
% Authored by T. Li, E-mail. tiancheng.li1985@gmail.com
% The variables have been explained in detail according to my paper: 
% T. Li, S. Sun, M. Bolic and J. M. Corchado. 
% Algorithm design for parallel implementation of the SMC-PHD filter, 
% Signal Processing, 2016, vol.119, pp. 115-127.
% Website: https://sites.google.com/site/tianchengli85/
%% ================== prediction of previous ===================
Dis = Y.Dis;
bearx = Y. bearing;
N = size(px,2);
Z = size(Dis,2);

wPre = pw*Ps;
xPre = px;
    Particle_State_Prediction
    y = xPre([1,3],:);
%% ================== add new target =====================
% L = ceil(Birth.w*Np);       % This serves the goal, or alternatively:
L= max(round(Birth.w*Np),200);% to make sure at least 100 particles used
sumL = sum(L);
        glik=zeros(N+sumL,Z);
if sumL >0
    for i= 1:numel(L)
        xPre=[xPre,GENRANDN(Birth.m(:,i),Birth.P,L(i))]; % one suggestion is to use a larger Birth.P?
        wPre=[wPre,sum(Birth.w(i))/L(i)*ones(1,L(i))]; 
    end
    for j=N+1:N+sumL
       y(:,j)=H*xPre(1:4,j); 
    end
end
%% ============== likelihood calculation of all =================
     %range
          Disy=sqrt((y(1,:)-sensor(1)).^2+(y(2,:)-sensor(2)).^2);
     %bearing     
          bearpxhat=atan2(y(2,:)-sensor(2),y(1,:)-sensor(1));
     for z=1:Z
          bearlik=normpdf(bearpxhat,bearx(z),r_noise);
          positionlik=normpdf(Disy,Dis(z),Wd);    
          glik(:,z)=bearlik.*positionlik;
     end
