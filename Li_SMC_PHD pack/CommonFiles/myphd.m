function [pw,wlik,Cz]=myphd(glik,Kz,Pd,wPre)    

% The core of Mahler's PHD updating 
%% input;
% glik:   The likelihood matrix, appearing as glik(j,z) 
%         where j index particle, and z index  measurements. 
%         glik(j,z) is corresponding to the g(z|j) in my paper;
% Kz:     The intensity of clutter, i.e. K_k(z)
% Pd:     The probability of detection, i.e. P_D in the paper
% wPre:   The weight of particles of the previous iteration, i.e. w_(k-1) 
%% output;
% pw:      Updated weight of particles, i.e. w_(k|k)
% wlik:  Weight-component matrix,appearing as wComp(j,z) 
%         where j index particle order, and z index the measurements. 
%         wComp(j,z) is w_k(z,j)in the paper
% Cz:     i.e. C_k(z) in the paper

% Writen by T. Li, 10/04/2013, e-mail. tiancheng.li1985@gmail.com

% The variables have been explained in detail according to my paper: 
% T. Li, S. Sun, M. Bolic and J. M. Corchado. 
% Algorithm design for parallel implementation of the SMC-PHD filter, 
% Signal Processing, 2016, vol.119, pp. 115-127.

% see also Multi-EAP: Extended EAP for multiple estimate extraction 
% for the SMC-PHD filter, submitted to IEEE TAES, 2015

% Preprint available @ Website: https://sites.google.com/site/tianchengli85/

%%
J=size(glik,1); % the number of particles
Z=size(glik,2); % the number of measurements
Cz=zeros(1,Z);  
    pw=zeros(1,J); 
    for z=1:Z
        Cz(z)=wPre*(glik(:,z).*Pd);% clik(j,z) i.e. c_k(z,j) in the paper;
    end
    wlik=zeros(J,Z);
    if numel(Pd)==1 %%% revised 25/2/2014, for varying detection probability
        for j=1:J
             wlik(j,:)=wPre(j)*glik(j,:)*Pd./(Kz+Cz); 
             pw(j)=(1-Pd)*wPre(j)+ sum(wlik(j,:));
        end
    else
        for j=1:J
             wlik(j,:)=wPre(j)*glik(j,:)*Pd(j)./(Kz+Cz); 
             pw(j)=(1-Pd(j))*wPre(j)+ sum(wlik(j,:));
        end
    end
 %%   old version for constant Pd
% J=size(glik,1); % the number of particles
% Z=size(glik,2); % the number of measurements
% clik=glik*Pd;   % clik(j,z) i.e. c_k(z,j) in the paper
%     Cz=zeros(1,Z);  
%     w=zeros(1,J); 
%     for z=1:Z
%         Cz(z)=wPre*clik(:,z);
%     end
%     wComp=zeros(J,Z);
%     for j=1:J
% %         for z=1:Z  %-- improved 24/11/2013---%
%             wComp(j,:)=wPre(j)*clik(j,:)./(Kz+Cz); 
% %         end
%          w(j)=(1-Pd)*wPre(j)+ sum(wComp(j,:));
%     end
    