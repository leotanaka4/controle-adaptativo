clc; clear; close all;

%% Definir Planta
kp = 0.1;
b = 3;
P = tf(kp*[1 2 1], [1 5 3 4]);

%% Definir Modelo
km = 1;
bm = 2;
M = tf(km*[1 3], [1 2 1]);

%% Polinômio do Observador
A0 = 1; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n,Lambda] = controle2DOF(P, M, A0);
