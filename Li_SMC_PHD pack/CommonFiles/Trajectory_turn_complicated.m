
function Targets=Trajectory_turn_complicated(F, G,Noise,Birth,Runtime,Window,sensor)

% nearly constant turn-rate (NCT) models with 10 targets born in 4 areas

% %     Birth.m = [-1500,0,250,0, 0; ...
% %         -250, 0, 1000, 0, 0;...
% %         250,0,750,0,0; ...
% %         1000,0,1500,0,0]'; 
x1=GENRANDN(Birth.m(:,1),Birth.P,1);
x2=GENRANDN(Birth.m(:,1),Birth.P,1);
x3=GENRANDN(Birth.m(:,2),Birth.P,1);
x4=GENRANDN(Birth.m(:,2),Birth.P,1);
x5=GENRANDN(Birth.m(:,3),Birth.P,1);
x6=GENRANDN(Birth.m(:,3),Birth.P,1);
x7=GENRANDN(Birth.m(:,3),Birth.P,1);
x8=GENRANDN(Birth.m(:,4),Birth.P,1);
x9=GENRANDN(Birth.m(:,4),Birth.P,1);
x10=GENRANDN(Birth.m(:,4),Birth.P,1);
% x1=[-1500,0,250,0, randn * Noise.turn]';
% x2=[-1500,0,250,0, randn * Noise.turn]';
% x3=[-250, 0, 1000, 0, randn * Noise.turn]';
% x4=[-250, 0, 1000, 0, randn * Noise.turn]';
% x5=[250,0,750,0, randn * Noise.turn]';
% x6=[250,0,750,0, randn * Noise.turn]';
% x7=[250,0,750,0, randn * Noise.turn]';
% x8=[1000,0,1500,0, randn * Noise.turn]';
% x9=[1000,0,1500,0, randn * Noise.turn]';
% x10=[1000,0,1500,0, randn * Noise.turn]';

p1=1;p2=p1;p3=p1;p4=p1;p5=p1;p6=p1;p7=p1;p8=p1;p9=p1;p10=p1;
X1=[];X2=[];X3=[];X4=[];X5=[];X6=[];X7=[];X8=[];X9=[];X10=[];
    for t=1:Runtime
        x=[];
        if t>=1 && t<=67 && p1
            if InsideWindow(x1,Window,sensor)
                x=[x,x1];
                X1=[X1,x1];
            else
                p1=0;
            end
%             w1 = w1 + randn*u;
            x1=F(x1) + [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 2===================
        if t>=10 && t<=58 && p2
            if InsideWindow(x2,Window,sensor)
                x=[x,x2];
                X2=[X2,x2];
            else
                p2=0;
            end
%             w2 = w2 + randn*u;
            x2=F(x2) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 3===================
        if t>=1 && t<=39 && p3
            if InsideWindow(x3,Window,sensor)
                x=[x,x3];
                X3=[X3,x3];
            else
                p3=0;
            end;
%             w3 = w3 + randn*u;
            x3=F(x3) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 4===================
        if t>=15 && t<=43 && p4
            if InsideWindow(x4,Window,sensor)
                x=[x,x4];
                X4=[X4,x4];
            else
                p4=0;
            end
%             w4 = w4 + randn*u; 
            x4=F(x4) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
            %================trajectory of target 5===================
        if t>=30 && t<=98 && p5
            if InsideWindow(x5,Window,sensor)
                x=[x,x5];
                X5=[X5,x5];
            else
                p5=0;
            end
%             w5 = w5 + randn*u;
            x5=F(x5) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 6===================

        if t>=20 && t<=80 && p6
            if InsideWindow(x6,Window,sensor)
                x=[x,x6];
                X6=[X6,x6];
            else
                p6=0;
            end
%             w6 = w6 + randn*u;
            x6=F(x6) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 7===================
        if t>=40 && t<=80 && p7
            if InsideWindow(x7,Window,sensor)
                x=[x,x7];
                X7=[X7,x7];
            else
                p7=0;
            end
%             w7 = w7 + randn*u;
            x7=F(x7) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 8===================
        if t>=39 && t<=80 &&p8
            if InsideWindow(x8,Window,sensor)
                x=[x,x8];
                X8=[X8,x8];
            else
                p8=0;
            end;
%             w8 = w8 + randn*u;
            x8=F(x8) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 9===================
        if t>=40 && t<=100 && p9
            if InsideWindow(x9,Window,sensor)
                x=[x,x9];
                X9=[X9,x9];
            else
                p9=0;
            end
%             w9 = w9 + randn*u; 
            x9=F(x9) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end
        %================trajectory of target 10===================
        if t>=60 && t<=100 && p10
            if InsideWindow(x10,Window,sensor)
                x=[x,x10];
                X10=[X10,x10];
            else
                p10=0;
            end
%             w10 = w10 + randn*u;
            x10=F(x10) +  [G*(Noise.vel.*randn(2,1)); randn * Noise.turn];
        end      
        Targets(t).x = x;
    end

filename= 'Trajectory_turn.mat';
save(filename,'Targets','X1','X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10');

function d=InsideWindow(xp,Window,sensor)
    d=0;
    dis=sqrt((xp(1)-sensor(1))^2+(xp(3)-sensor(2))^2);
    bearing=atan2(xp(3)-sensor(2),xp(1)-sensor(1));
if dis>=Window(1)&&dis<=Window(2)&&bearing>=Window(3)&&bearing<=Window(4)
    d=1;
end
