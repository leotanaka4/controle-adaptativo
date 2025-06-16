clc; close all; 
clear; 

%% ======= Escolher Ordem do Sistema =======
n = 2;  % ordem da planta e do modelo
reset = 1; % resetar a planta / usar a que já foi criada

%% ======= Inicialization =======
gamma = 10*eye(2*n);
tfinal = 300;    %Simulation interval
st = 0.01;      %Sample time to workspace
theta0 = zeros(2*n, 1);

%% ======= Reference signal parameters =======
DC = 1   %Constant
Aq = 2   %Sine wave amplitude
wq = 0.1*pi  %Frequency
As = 0   %Square wave amplitude
% As = 0.5
ws = pi  %Frequency


%% ======= Gerar Planta =======
if reset
    % Polos da planta (podem ser instáveis)
    plant_polos = randi([-10, 0], 1, n);           % polos arbitrários
    plant_den = poly(plant_polos);

    % Zeros reais negativos para numerador (mínima fase)
    plant_zeros = -randi([1, 5], 1, n-1);          % zeros reais negativos

    kp = 2 + rand();                               % ganho da planta
    plant_num = kp * poly(plant_zeros);            % numerador grau n-1

    P = tf(plant_num, plant_den);
    x0 = zeros(n, 1);
end
%% ======= Gerar Modelo SPR com polos e zeros alternados =======

% Passo 1: gerar polos reais negativos crescentes
pole_spacing = 2.0;
poles_M = sort(-1 * (1:pole_spacing:n*pole_spacing) );

% Passo 2: gerar zeros entre os polos (alternados)
zeros_M = zeros(1, n-1);
for i = 1:n-1
    zeros_M(i) = (poles_M(i) + poles_M(i+1)) / 2;  % entre dois polos
end

% Polinômios
den_M = poly(poles_M);
num_M = poly(zeros_M);

% Ganho positivo arbitrário
km = 1 + rand();
num_M = km * num_M;

M = tf(num_M, den_M);  % Modelo SPR garantido
xm0 = zeros(n,1);

%% Polinômio do Observador
A0 = 1; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0);
theta_star = [theta1; theta_n; theta2; theta_2n];

%% ======= Definir Filtro =======
% Grau da planta (denominador de P)
nf = n - 1;                 % grau do filtro

% Definir coeficientes do denominador do filtro (escolha arbitrária estável)
% Exemplo: filtro com polos em -2, -3, ..., arbitrário mas estável
filtro_polos = -2:-1:-(nf+1);                % nf polos reais negativos
den_filtro = poly(filtro_polos);            % coeficientes do denominador

% Forma canônica controlável
Af = zeros(nf); Af(1:nf-1, 2:nf) = eye(nf-1); Af(end, :) = -den_filtro(2:end);
Bf = zeros(nf,1); Bf(end) = 1;

%% ======= Simulação e Processamento =======
sim('MRAC_generalizado', tfinal);

% Parâmetros ideais replicados no tempo
theta_star_t = repmat(theta_star, 1, length(t));

% Normas dos parâmetros
norm_theta = vecnorm(theta');
norm_theta_star = vecnorm(theta_star_t);

% Entrada ideal
u_star = theta_star' * w';
u_star_t = repmat(u_star, 1, length(t));

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
plot(t, theta, 'LineWidth', 1.2); hold on;
plot(t, theta_star_t, '--', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('\theta');
title('\theta estimado vs \theta^* (ideal)');
legend(arrayfun(@(i)sprintf('\\theta_{%d}',i),1:2*n,'UniformOutput',false));
grid on;

% 3. Erro de acompanhamento
subplot(3,2,3);
plot(t, ea, 'm', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Erro');
title('Erro de acompanhamento: e_a = y_p - y_m');
grid on;

% 4. Norma de theta e theta*
subplot(3,2,4);
plot(t, norm_theta, 'b', 'LineWidth', 1.2); hold on;
plot(t, norm_theta_star, 'k--', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('||\theta||');
legend('||\theta||', '||\theta^*||');
title('Norma dos parâmetros');
grid on;

% 5. Controle u e u*
subplot(3,2,5);
plot(t, u, 'r', 'LineWidth', 1.2); hold on;
plot(t, u_star, 'b--', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('u(t)');
legend('u','u^*');
title('Controle aplicado vs ideal');
grid on;

% 6. Norma de theta tilde
subplot(3,2,6);
plot(t, norm_tilde, 'k', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('$$\|\tilde{\theta}\|$$','Interpreter','latex');
title('Erro de estimação: ||\theta - \theta^*||');
grid on;

