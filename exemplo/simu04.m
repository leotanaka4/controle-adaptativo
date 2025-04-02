%---------------------------------------------------------------------
%  COE-835  Controle adaptativo
%
%  Script para simular o exemplo 1
%                                                       Ramon R. Costa
%                                                       01/out/09, Rio
%---------------------------------------------------------------------
clear all;
clc;

%---------------------------------------------------------------------
disp('-------------------------------')
disp('Script para simular o exemplo 1')
disp(' ')
disp('Caso: Planta ............. n = 1')
disp('      Grau relativo ..... ns = 1')
disp('      Parâmetros ........ np = 1')
disp(' ')
disp('Algoritmo: Standard MRAC')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 10;    %Simulation interval
st = 0.05;      %Sample time to workspace

s = tf('s');    %trick!

%--------------------------------------------------------- Plant -----
ap = 2

P = 1/(s-ap)
P = ss(P);

%----------------------------------------------- Reference model -----
am = 1

M = 1/(s+am)
M = ss(M);

%--------------------------------------------- Initial condition -----
yp0  = 2
x0   = yp0;

ym0  = 0
xm0  = ym0;

%----------------------------------- Reference signal parameters -----
DC = 1   %Constant

As = 0   %Sine wave amplitude
ws = 10  %Frequency

%------------------------------------------------- Matching gain -----
thetas = -ap - am   %theta*

%----------------------------------------- Adaptation parameters -----
gamma = 10;      %Adaptation gains
theta0 = 0;      %Adaptation inicial condition

%---------------------------------------------- First simulation -----
[time,X] = sim('MRAC_111',tfinal);

T = t;     %Salva simulação
YP = yp;
YM = ym;
R = r;
E0 = e0;
Theta = theta;
U = u;

%---------------------------------------------------- First plot -----
figure(1)
clf
subplot(211)
hold on
grid on
f1 = plot(t,e0,'Color',[0 0.45 0.75],'Linew',0.5);
%style = ['Linew' num2str(0.5)]
%plot(t,e0,style);

Thetas = thetas*ones(size(T));

figure(2)
clf
subplot(211)
hold on
grid on
plot(t,Theta,'Color',[0 0.45 0.75],'Linew',0.5);
plot(t,Thetas,'Color',[0.85 0.35 0.1],'Linew',0.5);
%f2 = plot(t,Theta,'b',t,Thetas,'r','Linew',0.5);
%f2.Color

figure(3)
clf
subplot(211)
hold on
grid on
plot(t,yp,'Color',[0 0.45 0.75],'Linew',0.5);
plot(t,ym,'Color',[0.85 0.35 0.1],'Linew',0.5);

%---------------------------------------- Continuous simulation -----
Ns = 5       %Número de simulações
Dt = 10;     %Incremento no tempo de simulação
Dg = 0;      %Incremento de gamma
Da = 1;      %Incremento no parâmetro ap

for i=2:Ns
    Xfinal = X(size(X,1),:);        %Estado final da última simulação
    X0 = Xfinal;             %Estado inicial para a próxima simulação

    t0 = tfinal;      
    tfinal = t0 + Dt;
    intervalo = [t0 tfinal];     %Intervalo da próxima simulação
    
    gamma = gamma + Dg;          %Alteração de parâmetros
    ap = ap + Da
    P = ss(1/(s-ap));
    thetas = -ap - am   %theta*

    %opcoes = simset('StartTime',num2str(t0),'StopTime',num2str(tfinal),'InitialState',num2str(X0));
    opcoes = simset('InitialState',X0);
    
    %[time,X] = sim('MRAC_111',intervalo,opcoes); %Continuação!
    [time,X] = sim('MRAC_111',intervalo,opcoes); %Continuação!
    %[time,X] = sim('MRAC_111','StartTime',t0,'StopTime',tfinal,'InitialState',X0)

    T = [T ; t];     %Save results
    YP = [YP ; yp];
    YM = [YM ; ym];
    R = [R ; r];
    E0 = [E0 ; e0];
    Theta = [Theta ; theta];
    U = [U ; u];

    figure(1)
    plot(t,e0,'Color',[0 0.45 0.75],'Linew',0.5);

    Thetas = thetas*ones(size(theta));
    
    figure(2)
    plot(t,theta,'Color',[0 0.45 0.75],'Linew',0.5);
    plot(t,Thetas,'Color',[0.85 0.35 0.1],'Linew',0.5);

    figure(3)
    plot(t,yp,'Color',[0 0.45 0.75],'Linew',0.5);
    plot(t,ym,'Color',[0.85 0.35 0.1],'Linew',0.5);
end

Thetas = thetas*ones(size(T));

%----------------------------------------------- Print eps plots -----
figure(1)
subplot(211)
title({'$e_0$'},'FontSize',10,'Interpreter','latex')
legend({'$e_0$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig04a.eps

figure(2)
subplot(211)
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
legend({'$\theta$','$\theta^*$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig04b.eps

figure(3)
subplot(211)
title({'$y_m, y_p$'},'FontSize',10,'Interpreter','latex')
legend({'$y_m$','$y_p$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig04c.eps

figure(4)
clf
subplot(211)
plot(T,Theta,T,Thetas,'Linew',0.5);
grid on; 
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
legend({'$\theta$','$\theta^*$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig04d.eps

%------------------------------------------------- Display plots -----
figure(5)
clf

subplot(221)
plot(T,E0,'Linew',0.5);
grid on; 
title({'$e_0$'},'FontSize',10,'Interpreter','latex')
legend({'$e_0$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')

subplot(222)
plot(T,Theta,T,Thetas,'r','Linew',0.5);
grid on; 
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
legend({'$\theta$','$\theta^*$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')

subplot(223)
plot(T,YP,T,R,T,YM,'Linew',0.5);
grid on; 
title({'$y_m, y_p$'},'FontSize',10,'Interpreter','latex')
legend({'$y_m$','$y_p$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')

subplot(224)
plot(T,U,'Linew',0.5);grid;
grid on; 
title({'$u$'},'FontSize',10,'Interpreter','latex')
legend({'$u$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')

%--------------------------------------- Impressão dos diagramas -----
% open_system('MRAC_111');
% print -depsc2 -sMRAC_111 MRAC-111.eps
% 
% open_system('MRAC_111/Plant');
% print -depsc2 -sMRAC_111/Plant plant.eps
% 
% open_system('MRAC_111/Reference model');
% print -depsc2 '-sMRAC_111/Reference model' reference-model.eps
% 
% open_system('MRAC_111/Adaptation');
% print -depsc2 -sMRAC_111/Adaptation adaptation.eps
% 
% open_system('MRAC_111/Reference signal');
% print -depsc2 '-sMRAC_111/Reference signal' reference-signal.eps
% 
% close_system('MRAC_111');
%---------------------------------------------------------------------


