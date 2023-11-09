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


%% simulate a network
%clear all
%close all
tspan = [0 100];
%flux =ones(10,1); %randomly generate fluxes
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
%y0 = ones(6,1); %randomly generate starting metabolites
y0 = ones(6,1);

[t,y] = ode45(@(t,y) LIFEsimRCT(y,flux), tspan, y0);

%% Plot metabolite trajectories
figure(1)
plot(t,y,'LineWidth',2)
title('Metabolite values over time','FontSize',20)
ylabel('X','FontSize',20)
xlabel('time','FontSize',20)

%% derivative function
function derivative = LIFEsimRCT(x,flux)
    S = [1 0 0 -x(1) 0 0 0 0 0 0;
        0 1 0 0 -x(2) 0 0 0 0 0;
        0 0 1 0 0 -x(3) 0 0 0 0;
        0 0 0 x(1) x(2) x(3) -x(4) -x(4) 0 0;
        0 0 0 0 0 0 x(4) 0 -x(5) 0;
        0 0 0 0 0 0 0 x(4) x(5) -x(6)];

    derivative = S*flux;
end
