clc; clear; close all;

%% Definir Planta
kp = 2;
b = 3;
P = tf([kp b*kp], [1 4 3]);

%% Definir Modelo
km = 1;
bm = 2;
M = tf([km bm*km], [1 5 6]);

%% Polinômio do Observador
A0 = 1; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0);
