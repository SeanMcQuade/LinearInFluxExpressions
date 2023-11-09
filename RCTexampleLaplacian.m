% display latex

X=sym('x', [1,6]);
A=sym('a', [1,4]);
S = [1 0 0 -X(1) 0 0 0 0 0 0;
    0 1 0 0 -X(2) 0 0 0 0 0;
    0 0 1 0 0 -X(3) 0 0 0 0;
    0 0 0 X(1) X(2) X(3) -X(4) -X(4) 0 0;
    0 0 0 0 0 0 X(4) 0 -X(5) 0;
    0 0 0 0 0 0 0 X(4) X(5) -X(6)];
N_S = null(S)


sum_Null = A(1)*N_S(:,1) + A(2)*N_S(:,2) + A(3)*N_S(:,3) + A(4)*N_S(:,4)
latex(S)
latex(N_S(:,1))
latex(N_S(:,2)) 
latex(N_S(:,3)) 
latex(N_S(:,4)) 

latex(sum_Null)



flux =[0.7017
    0.6479
    0.7520
    0.9196
    0.7402
    0.2922
    0.8157
    0.4870
    0.5868
    0.8255];
%% calculate equilibria
%J*x=-phi -> 
%J_times_x = J*x


%phi = transpose([f_1, f_2, f_3, 0, 0, 0]);
% L = ['-f_4', '0', '0', '0', '0', '0';
%      '0', '-f_5', '0', '0', '0', '0';
%      0, 0, '-f_6', 0, 0, 0;
%      'f_4','f_5','f_6','-(f_7+f_8)',0,0;
%      0, 0, 0, 'f_8', '-f_9', 0;
%      0, 0, 0, 'f_7', 'f_9', '-f_{10}'];

% L = ['-f_4', '0', '0', '0', '0', '0';
%      '0', '-f_5', '0', '0', '0', '0';
%      '0', '0', '-f_6', '0', '0', '0';
%      'f_4','f_5','f_6','-(f_7+f_8)','0','0';
%      '0', '0', '0', 'f_8', '-f_9', '0';
%      '0', '0', '0', 'f_7', 'f_9', '-f_{10}'];

f=sym('f', [10,1]);
x=sym('x', [6,1]);
phi = transpose([f(1), f(2), f(3), 0, 0, 0]);
latex(phi)
L = [-f(4), 0, 0, 0, 0, 0;
     0, -f(5), 0, 0, 0, 0;
     0, 0, -f(6), 0, 0, 0;
     f(4),f(5),f(6),-(f(7)+f(8)),0,0;
     0, 0, 0, f(8), -f(9), 0;
     0, 0, 0, f(7), f(9), -f(10)];



%% simulate a network
%clear all
%close all
tspan = [0 50];
y0 = 4*rand(6,1); %randomly generate starting metabolites
%y0 = [2;2;1;1;1;2];
%flux = ones(10,1); %randomly generate fluxes


[t,y] = ode45(@(t,y) LIFEsimRCT(y,flux), tspan, y0);



%% see derivative
flux =[0.7017
    0.6479
    0.7520
    0.9196
    0.7402
    0.2922
    0.8157
    0.4870
    0.5868
    0.8255];
phi = transpose([flux(1), flux(2), flux(3), 0, 0, 0]);
L = [-flux(4), 0, 0, 0, 0, 0;
     0, -flux(5), 0, 0, 0, 0;
     0, 0, -flux(6), 0, 0, 0;
     flux(4),flux(5),flux(6),-(flux(7)+flux(8)),0,0;
     0, 0, 0, flux(8), -flux(9), 0;
     0, 0, 0, flux(7), flux(9), -flux(10)];

initial_derivative = L*y0+phi;
final_derivative = L*transpose(y(end,:))+phi;
equilibria = -inv(L)*phi;

%% Plot metabolite trajectories
figure(2)
plot(t,y,'LineWidth',2)
title('Reverse Cholesterol Transport','FontSize',20)
ylabel('Metabolite levels X','FontSize',20)
xlabel('time (hours)','FontSize',20)
xlim(tspan)
ylim([0.5 4])
hold on
plot(linspace(0,50),equilibria(1)*ones(100,1), 'Linestyle','--','LineWidth',2)
hold on
plot(linspace(0,50),equilibria(2)*ones(100,1), 'LineStyle','--','LineWidth',2)
hold on
plot(linspace(0,50),equilibria(3)*ones(100,1), 'LineStyle','--','LineWidth',2)
hold on
plot(linspace(0,50),equilibria(4)*ones(100,1), 'LineStyle','--','LineWidth',2)
hold on
plot(linspace(0,50),equilibria(5)*ones(100,1), 'LineStyle','--','LineWidth',2)
hold on
plot(linspace(0,50),equilibria(6)*ones(100,1), 'LineStyle','--','LineWidth',2)
% yyaxis right
% ylim([0 4])
 %yticks([0.7630, 0.8753, 1.3389, 1.6133, 2.5459, 2.5736])
%% derivative function
function derivative = LIFEsimRCT(y,flux)
% LIFE dyn
%     S = [1 0 0 -y(1) 0 0 0 0 0 0;
%         0 1 0 0 -y(2) 0 0 0 0 0;
%         0 0 1 0 0 -y(3) 0 0 0 0;
%         0 0 0 y(1) y(2) y(3) -y(4) -y(4) 0 0;
%         0 0 0 0 0 0 y(4) 0 -y(5) 0;
%         0 0 0 0 0 0 0 y(4) y(5) -y(6)];
% 
%     derivative = S*flux;

%Laplacian dyn
phi = transpose([flux(1), flux(2), flux(3), 0, 0, 0]);
L = [-flux(4), 0, 0, 0, 0, 0;
     0, -flux(5), 0, 0, 0, 0;
     0, 0, -flux(6), 0, 0, 0;
     flux(4),flux(5),flux(6),-(flux(7)+flux(8)),0,0;
     0, 0, 0, flux(8), -flux(9), 0;
     0, 0, 0, flux(7), flux(9), -flux(10)];
derivative = L*y+phi;
end