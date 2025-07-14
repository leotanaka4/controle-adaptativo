clc; close all; 
clear all;
%% ======= Escolher Ordem do Sistema =======
n = 3;  % ordem da planta e do modelo

%% ======= Inicialization =======
y0 = 0;
ym0 = 0;
gamma = 20;
R0 = 20*eye(2*n);
tfinal = 150;    %Simulation interval
st = 0.01;      %Sample time to workspace
theta0 = zeros(2*n, 1);

%% ======= Reference signal parameters =======
DC = 3;   %Constant
S1 = 0;
S2 = 0;
S3 = 0;
SQ = 1;
% Define Laplace variable
s = tf('s');
kp= 1/3;
P = kp*((s + 2)^2) / (s^3);
P = ss(P);
x0 = P.c \ y0;
M = 1 / (s + 1);
M = ss(M);
xm0 = M.c \ ym0;
Lambda = 1 / ((s + 0.5)*(s + 1)*(s + 1.5));
l0 = 2;

%% Polinômio do Observador Vetorizado
A0 = 1; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n, den_filtro] = controle2DOF(P, M, A0);
theta_star = [theta1(:); theta_n; theta2(:); theta_2n];
theta0 = 0.9*theta_star;
%% ======= Definir Filtro =======
% Forma canônica controlável
nf = n-1;
Af = zeros(nf); Af(1:nf-1, 2:nf) = eye(nf-1); Af(end, :) = -fliplr(den_filtro(2:end));
Bf = zeros(nf,1); Bf(end) = 1;
%theta0 = theta_star;

sim('LS_MRAC', tfinal);

% Parâmetros ideais replicados no tempo
theta_star_t = repmat(theta_star, 1, length(t));

% Normas dos parâmetros
norm_theta = vecnorm(theta');
norm_theta_star = vecnorm(theta_star_t);

% Erro de parâmetros e sua norma
theta_til = theta' - theta_star_t;
norm_tilde = vecnorm(theta_til);

%% === Gráficos ===
figure('Name','Resultados da Simulação MRAC','Units','normalized','Position',[0.05 0.05 0.9 0.85]);

% 1. Referência, ym, yp
subplot(3,2,1);
plot(t, r, 'k', 'LineWidth', 1.2); hold on;
plot(t, ym, 'b--', 'LineWidth', 1.2);
plot(t, yp, 'r', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Saída');
legend('r','y_m','y_p'); title('Saídas: Referência, Modelo, Planta');
grid on;

% 2. theta vs theta*
subplot(3,2,2);
colors = lines(2*n);
hold on;
for i = 1:2*n
    plot(t, theta(:,i), 'Color', colors(i,:), 'LineWidth', 1.2);
    plot(t, theta_star_t(i,:), '--', 'Color', colors(i,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
xlabel('Tempo [s]');
ylabel('\theta');
title('\theta estimado vs \theta^* (ideal)');
legend(arrayfun(@(i)sprintf('\\theta_{%d}', i), 1:4*n, 'UniformOutput', false));
grid on;

% 3. Erro de acompanhamento
subplot(3,2,3);
plot(t, ea, 'm', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Erro');
title('Erro de acompanhamento: e_a = y_p - y_m');
grid on;

% 4. Norma de theta e theta*
subplot(3,2,5);
plot(t, norm_theta, 'b', 'LineWidth', 1.2); hold on;
plot(t, norm_theta_star, 'k--', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('||\theta||');
legend('||\theta||', '||\theta^*||');
title('Norma dos parâmetros');
grid on;

% 5. Controle u e u*
subplot(3,2,4);
plot(t, sum(theta_til.*w',1) + sum(dtheta'.*w',1), 'LineWidth', 1.2); hold on;
plot(t, sum(theta_til.*w',1), 'LineWidth', 1.2);
plot(t, sum(dtheta'.*w',1), 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('u(t)');
legend('$\tilde{u}$', '$\tilde{\theta}^T \omega$', '$\dot{\theta}^T \omega$', 'Interpreter', 'latex');

title('Controle aplicado vs ideal');
grid on;

% 6. Norma de theta tilde
subplot(3,2,6);
plot(t, norm_tilde, 'k', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('$$\|\tilde{\theta}\|$$','Interpreter','latex');
title('Erro de estimação: ||\theta - \theta^*||');
grid on;
