
function main_MEAP

% the SMC-PHD filter, based on nearly constant turn-rate (NCT) model.

% The multi-estimate extraction is based on the so-called MEAP estimator
% which is much more accurate and fater than the commonly used clustering

% For the MEAP method for multi-estimate extraction,see:
% T. Li, J. M. Corchado, S. Sun and H. Fan, Multi-EAP: Extended EAP for 
% multiple estimate extraction  for the SMC-PHD filter, Chinese Journal
% of Aeronautics, to apprear.
% See also: T. Li, S. Sun, M. Bolic and J. M. Corchado. 
% Algorithm design for parallel implementation of the SMC-PHD filter, 
% Signal Processing, 2016, vol.119, pp. 115-127.

% Contact: tiancheng.li1985@gmail.com;
% Website: https://sites.google.com/site/tianchengli85/matlab-codes/phd-filter

addpath ./CommonFiles

%% parameter setting
%%
close all
clear

T=1;%sampling period
F = @(a)[1 sin(a(5,:)*T)/a(5,:)      0    -((1-cos(a(5,:)*T))/a(5,:)) 0
         0 cos(a(5,:)*T)             0    -sin(a(5,:)*T)              0
         0 (1-cos(a(5,:)*T))/a(5,:)  1    sin(a(5,:)*T)/a(5,:)        0
         0 sin(a(5,:)*T)             0    cos(a(5,:)*T)               0
         0       0                   0    0                           1]*a;
    Noise.vel = [15; 15];
    Noise.turn= pi/180;

    G= [T^2/2 0;T 0;0 T^2/2;0 T];
    H = [1 0 0 0;0 0 1 0];
    Wd = 5; % Wd = 10;
    sensor=[0;0]; 
    r_noise= pi/180; % r_noise=2*pi/180 

% %  Parametes setting  
    Model.pos= H;
    Model.Pd=0.95;
    Ps=0.99;
    Rate_Clutter = 10; % clutter avarage number  
    Kz=Rate_Clutter/2000/pi;
%  miss distance  by OPFA
    w_Ristic=0.6;
    c=500; %  it is better to use a cutoff parameter related to the scenario
    p=2;
    Np=1000; %  number of particle allocated to each expected target
    Nmin=600;      
%%% based on fixed position, Vo's model
    Birth.m = [-1500,0,250,0, 0; ...
        -250, 0, 1000, 0, 0;...
        250,0,750,0,0; ...
        1000,0,1500,0,0]';             %Birth denote the apprear density
    Birth.P = diag([50 50 50 50 6*(pi/180)])^2;
    Birth.w = [0.02 0.02 0.03 0.03];
    %==============target trajectory===============
    Runtime= 100;
    Window =[0, 2000, 0, pi];
    Num_Targets(1)= 0;
    Targets=Trajectory_turn_complicated(F, G,Noise,Birth,Runtime,Window,sensor);
    load Trajectory_turn
    plotTrajectory

%%    %==================================main procedure==================
repeatTimes =2;  % 
sumd_pf=[];
sumd_MEAP_Nk=[];
sumd_MEAP_adaptive=[];
sumd_Zhao=[];
sumd_Ristic=[];
sumTime_pf=[];
sumTime_MEAP_Nk=[];
sumTime_MEAP_adaptive=[];
sumTime_Zhao=[];
sumTime_Ristic=[];
sumTNo_PFs=[];
sumTNo_Ristic=[];
sumTno=[];
sumTNo_MEAP_adaptive=[];
for repe=1:repeatTimes
    % initialization
        L = max(round(Birth.w*Np),Nmin);
        sumL = sum(L);
        xPre1=GENRANDN(Birth.m(:,1),Birth.P,L(1));
        xPre2=GENRANDN(Birth.m(:,2),Birth.P,L(2));
        xPre3=GENRANDN(Birth.m(:,3),Birth.P,L(3));
        xPre4=GENRANDN(Birth.m(:,4),Birth.P,L(4));  
        px  = [xPre1,xPre2,xPre3,xPre4];
        pw =sum(Birth.w)/sumL*ones(1,sumL);
%         N=sumL;
       TNo_PFs=[];
       TNo_Ristict=[];
       TNo_MEAP_adaptive=[];
       Time_pf=[];
       Time_Zhao=[];
       Time_Ristic=[];
       Time_MEAP_Nk=[];
       Time_MEAP_adaptive=[];
    d_pf=[]; 
    d_Ristic=[]; 
    d_Zhao=[]; 
    d_MEAP_Nk=[]; 
    d_MEAP_adaptive=[]; 
    for t=1:Runtime
        x=Targets(t).x;
        Num_Targets(t)=size(x,2);  
         clutterNo = poissrnd(Rate_Clutter);
        Y.Dis = rand(1,clutterNo)*2000;
            Y.bearing = rand(1,clutterNo)*pi;  
        if ~isempty(x)
            distance= H*x(1:4,:);
            Jp= find(rand(size(distance,2),1)<compute_pD(Model,x));
            Y.Dis = [Y.Dis sqrt((distance(1,Jp)-sensor(1)).^2+(distance(2,Jp)-sensor(2)).^2)+ Wd*randn(1,numel(Jp))];
            Y.bearing = [Y.bearing atan2(distance(2,Jp)-sensor(2),distance(1,Jp)-sensor(1))+ r_noise*randn(1,numel(Jp))];  
        end   
        observation(t).Y = Y;
        
         if ~isempty(Y)  
%%%% ============== comperison===============%    difference one  0 or 1
%             general
            [pxhat,what,wPre] = CompWeight_distance_bearing(Y,pw,px,H,F,G,Noise,sensor,r_noise,Birth,Np,Wd,Ps);
%             w=phd(what,Kz,pD,wPre);
        pD=compute_pD(Model,pxhat);
            [pw,wlik,Cz]=myphd(what,Kz,pD,wPre); 
        %% resampling
        sumw=sum(pw);
        Tno=round(sumw);
        TNo_PFs(t)=Tno;
        N=max(round(Np*sumw),Nmin);   %  N=max(Np*Tno,Nmin);
                outIndex=resampling(pw,N); 
                pw=sumw/N*ones(1,N);
                px=pxhat(:,outIndex(:))+ [G*(repmat(Noise.vel,1,N).*randn(2,N));zeros(1,N)]; 
                % for the roughening noises, please refer to the following work
                % T. Li, et al.Roughening methods to prevent sample impoverishment in the particle PHD filter, 
                % FUSION 2013, Istanbul Turkey, 9-12 July 2013.

%%  ///// OUTPUT   ////////%%%  
%% Clustering method
tic
%  Standard k-means Clustering
            xout=cluster(px([1,3],:),Tno); 
    Time_pf(t)=toc;
%   xoutArr(t).xout=xout;
%% Ristic's method
tic
Tno2=0;
L=size(Y.Dis,2); 
        Wz=sum(wlik);
        Wz(L+1)=wPre*(1-pD);% sum((1-pD).*wPre);
            xout2=[];
            for z=1:L+1
                if Wz(z)>=w_Ristic
                    Tno2=Tno2+1;
                    if z==L+1
                        xout2=[xout2,pxhat([1,3],:)*(wPre'.*(1-pD))]; 
                    else
                        xout2=[xout2,pxhat([1,3],:)*wlik(:,z)];   
                    end
                end
            end
Time_Ristic(t)=toc;            
             TNo_Ristict(t)=Tno2; 
%% peak-extraction of Lingling Zhao's method
tic
L=size(Y.Dis,2);
Wz=sum(wlik);              
        Wz(L+1)=wPre*(1-pD);% sum();
        [Dtemp,wSub]=sort(-Wz);
%             xout3=zeros(2,Tno);
xout3=[];
            for jtemp=1:min(Tno,L+1)
                if wSub(jtemp)==L+1
                   wsub=wPre'.*(1-pD)/Wz(L+1);
                   xout3(:,jtemp)=pxhat([1,3],:)*wsub;   
                else
                    wsub=wlik(:,wSub(jtemp))/Wz(wSub(jtemp));
                    xout3(:,jtemp)=pxhat([1,3],:)*wsub; 
                end
            end
Time_Zhao(t)=toc;
%%   my method with the assumed number of estimates !
tic
L = size(Y.Dis,2);
Wz=sum(wlik);
            [Dtemp,wSub]=sort(-Wz);
            [Zmtx, indx] = max(wlik,[],2);
            xout4=zeros(2,min(Tno,L));
            for j=1:min(Tno,L)
                  Jz = find(indx==wSub(j));  % The nearest neighbor data association 
%                 Jz2 = find(what(:,wSub(j))>=0.058);% The 1/2/3 sigma near data association, 
%                  % where 0.058 = normpdf(1,0,1)^2; 0.0029 = normpdf(2,0,1)^2; 1.96e-5 = normpdf(3,0,1)^2;
%                 Jz = union(Jz,Jz2); %% Near and Nearest Neighbor data association.
%                                       % to note, Jz can choose just Jz1 or
%                                       % Jz2 for simplicity and more
%                                       % better estimation even!
                    wsubsum=wPre(Jz)*what(Jz,wSub(j));
                    wsub=wPre(Jz)'.*what(Jz,wSub(j))/wsubsum;
                  xout4(:,j)=pxhat([1,3],Jz)*wsub;
            end
Time_MEAP_Nk(t) = toc;           
%%  my method with adaptive number of estimates!!
tic
Tno1  = 0;
L = size(Y.Dis,2);
Wz=sum(wlik);
            [Zmtx, indx] = max(wlik,[],2);
                xout1=[];
                for z=1:L
                    if Wz(z)>=w_Ristic
                        Tno1=Tno1+1;
                        Jz= find(indx==z);
                        wsubsum=wPre(Jz)*what(Jz,z);
                        wsub=wPre(Jz)'.*what(Jz,z)/wsubsum;
                        xout1 = [xout1, pxhat([1,3],Jz)*wsub];
                    end
                end            
Time_MEAP_adaptive(t) = toc; 
TNo_MEAP_adaptive (t)  = Tno1;
%%         
            if ~isempty(x)
               d_pf(t)=ospa_dist(xout,x([1,3],:),c,p); 
               d_Ristic(t)=ospa_dist(xout2,x([1,3],:),c,p); 
               d_Zhao(t)=ospa_dist(xout3,x([1,3],:),c,p); 
               d_MEAP_Nk(t)=ospa_dist(xout4,x([1,3],:),c,p);
               d_MEAP_adaptive(t)=ospa_dist(xout1,x([1,3],:),c,p);
            else
               d_pf(t)=ospa_dist(xout,[],c,p); 
               d_Ristic(t)=ospa_dist(xout2,[],c,p); 
               d_Zhao(t)=ospa_dist(xout3,[],c,p);
               d_MEAP_Nk(t)=ospa_dist(xout4,[],c,p);
               d_MEAP_adaptive(t)=ospa_dist(xout1,[],c,p);
            end
        end
    end
    sumd_pf=[sumd_pf;d_pf];
    sumd_Ristic=[sumd_Ristic;d_Ristic];
    sumd_Zhao=[sumd_Zhao;d_Zhao];
    sumd_MEAP_Nk=[sumd_MEAP_Nk;d_MEAP_Nk];
    sumd_MEAP_adaptive=[sumd_MEAP_adaptive;d_MEAP_adaptive];
    sumTime_pf=[sumTime_pf;Time_pf];
    sumTime_Ristic=[sumTime_Ristic;Time_Ristic];
    sumTime_Zhao=[sumTime_Zhao;Time_Zhao];
    sumTime_MEAP_Nk=[sumTime_MEAP_Nk;Time_MEAP_Nk];
    sumTime_MEAP_adaptive=[sumTime_MEAP_adaptive;Time_MEAP_adaptive];
sumTNo_PFs=[sumTNo_PFs;TNo_PFs];
sumTNo_Ristic=[sumTNo_Ristic;TNo_Ristict];
sumTNo_MEAP_adaptive=[sumTNo_MEAP_adaptive;TNo_MEAP_adaptive];
sumTno=[sumTno;Num_Targets];

end
meand_pf=mean(sumd_pf);
meand_Ristic=mean(sumd_Ristic);
meand_Zhao=mean(sumd_Zhao);
meand_MEAP_Nk=mean(sumd_MEAP_Nk);
meand_MEAP_adaptive=mean(sumd_MEAP_adaptive);
meanTime_pf=mean(sumTime_pf);
meanTime_Ristic=mean(sumTime_Ristic);
meanTime_Zhao=mean(sumTime_Zhao);
meanTime_MEAP_Nk=mean(sumTime_MEAP_Nk);
meanTime_MEAP_adaptive=mean(sumTime_MEAP_adaptive);
meanTNo_PFs=mean(sumTNo_PFs);
meanTNo_Ristic=mean(sumTNo_Ristic);
meanTNo_MEAP_adaptive=mean(sumTNo_MEAP_adaptive);
meanTno=mean(sumTno);

save MEAP(test)

    Errorpecent_Ristic=(mean(meand_pf)-mean(meand_Ristic))/mean(meand_pf);
    Errorpecent_Zhao=(mean(meand_pf)-mean(meand_Zhao))/mean(meand_pf);
    Errorpecent_MEAP_Nk=(mean(meand_pf)-mean(meand_MEAP_Nk))/mean(meand_pf);
    Errorpecent_MEAP_adaptive=(mean(meand_pf)-mean(meand_MEAP_adaptive))/mean(meand_pf);

    Timepecent_Ristic=(mean(meanTime_pf)-mean(meanTime_Ristic))/mean(meanTime_pf);
    Timepecent_Zhao=(mean(meanTime_pf)-mean(meanTime_Zhao))/mean(meanTime_pf);
    Timepecent_MEAP_Nk=(mean(meanTime_pf)-mean(meanTime_MEAP_Nk))/mean(meanTime_pf);
    Timepecent_MEAP_adaptive=(mean(meanTime_pf)-mean(meanTime_MEAP_adaptive))/mean(meanTime_pf);

%%  ======== figures plot===========
%  mean OSPA plot
% close all
figure
subplot(3,1,1)
plot(meand_pf,'k-o'), 
hold on, plot(meand_Ristic,'g-*'),hold on, plot(meand_Zhao,'r-+'),
hold on, plot(meand_MEAP_adaptive,'m-x'),hold on, plot(meand_MEAP_Nk,'b.-'),
set(gca,'FontSize',12); set(gcf,'Color','White'); 
%h1=xlabel('Step');
h2=ylabel('Mean OSPA');
% h3=legend('{\itk}-means clustering','Ristic\primes method','Zhao\primes method','MEAP','location','best');%,'Orientation','horizontal'
set([h2],'FontSize',14,'FontName','Times new Roman');%set(h3,'FontSize',12,'FontName','Times new Roman');
% LEGEND BOXON ;
subplot(3,1,2)
semilogy(meanTime_pf,'k-o'), 
hold on, plot(meanTime_Ristic,'g-*'),hold on, plot(meanTime_Zhao,'r-+'),
hold on, plot(meanTime_MEAP_adaptive,'m-x'),hold on, plot(meanTime_MEAP_Nk,'b.-'),
set(gca,'FontSize',12); set(gcf,'Color','White'); 
% h1=xlabel('Step'); 
h2=ylabel('Computing time/s');
h3=legend('{\itk}-means clustering','Ristic\primes method','Zhao\primes method','MEAP I','MEAP II','location','best');
set([h2],'FontSize',14,'FontName','Times new Roman');set(h3,'FontSize',12,'FontName','Times new Roman');
% LEGEND BOXON ;
subplot(3,1,3)
plot(meanTno,'k-'),hold on, plot(meanTNo_Ristic,'g-*'),plot(meanTNo_MEAP_adaptive,'m-x'), hold on,plot(meanTNo_PFs,'r.-') 
set(gca,'FontSize',12); set(gcf,'Color','White'); 
h1=xlabel('Step'); h2=ylabel('Number of targets');
h3=legend('True number of targets','Ristic\primes method','MEAP I','Zhao\primes method and MEAP II','location','best');
set([h1,h2],'FontSize',14,'FontName','Times new Roman');set(h3,'FontSize',12,'FontName','Times new Roman');


LEGEND BOXON ;
figure(7)
subplot(2,1,1)
h1=xlabel('Step'); h2=ylabel('Distance');
hold on
    plot(1:size(X1,2),sqrt(X1(1,:).^2+X1(3,: ).^2),'b-','LineWidth',2),
    plot(10:9+size(X2,2),sqrt(X2(1,:).^2+X2(3,: ).^2),'b-','LineWidth',2),
    plot(1:size(X3,2),sqrt(X3(1,:).^2+X3(3,: ).^2),'b-','LineWidth',2),
    plot(15:14+size(X4,2),sqrt(X4(1,:).^2+X4(3,: ).^2),'b-','LineWidth',2)
    plot(30:29+size(X5,2),sqrt(X5(1,:).^2+X5(3,: ).^2),'b-','LineWidth',2)
    plot(20:19+size(X6,2),sqrt(X6(1,:).^2+X6(3,: ).^2),'b-','LineWidth',2),
    plot(40:39+size(X7,2),sqrt(X7(1,:).^2+X7(3,: ).^2),'b-','LineWidth',2)
    plot(39:38+size(X8,2),sqrt(X8(1,:).^2+X8(3,: ).^2),'b-','LineWidth',2)
    plot(40:39+size(X9,2),sqrt(X9(1,:).^2+X9(3,: ).^2),'b-','LineWidth',2)
    plot(60:59+size(X10,2),sqrt(X10(1,:).^2+X10(3,: ).^2),'b-','LineWidth',2)
subplot(2,1,2)
h3=xlabel('Step'); h4=ylabel('Bearing');
hold on
    plot(1:size(X1,2),atan2(X1(3,:),X1(1,:)),'b-','LineWidth',2),
    plot(10:9+size(X2,2),atan2(X2(3,:),X2(1,:)),'b-','LineWidth',2),
    plot(1:size(X3,2),atan2(X3(3,:),X3(1,:)),'b-','LineWidth',2),
    plot(15:14+size(X4,2),atan2(X4(3,:),X4(1,:)),'b-','LineWidth',2)
    plot(30:29+size(X5,2),atan2(X5(3,:),X5(1,:)),'b-','LineWidth',2)
    plot(20:19+size(X6,2),atan2(X6(3,:),X6(1,:)),'b-','LineWidth',2),
    plot(40:39+size(X7,2),atan2(X7(3,:),X7(1,:)),'b-','LineWidth',2)
    plot(39:38+size(X8,2),atan2(X8(3,:),X8(1,:)),'b-','LineWidth',2)
    plot(40:39+size(X9,2),atan2(X9(3,:),X9(1,:)),'b-','LineWidth',2)
    plot(60:59+size(X10,2),atan2(X10(3,:),X10(1,:)),'b-','LineWidth',2)
 axis([0,Runtime,0,3.5])
legend BOXOFF ;
set(gcf,'Color','White'); 
set([h1,h2,h3,h4],'FontSize',14,'FontName','Times new Roman');
%  observation plot
for t=1:Runtime
    Y =observation(t).Y ;
    if ~isempty(Y)
        figure(7)
        subplot(2,1,1)
        hold on
        h1=plot(t,Y.Dis,'o');set(h1,'Color','Black')
        subplot(2,1,2)
        hold on
        h2=plot(t,Y.bearing,'o');set(h2,'Color','Black')
    end
end