%% graus da planta
n = 5; %editar aqui - grau do denominador
m = 3; %editar aqui - grau do numerador
n_star = n-m; %editar aqui - g1rau relativo


%% inicializando sinais das plantas
u = 0;
r = 0;
y = 0;
e0 = 0;
ea - 0;

%% inicializando parametros das plantas
kp = 0;
km = 0;
bp = ones(1,m+1); 
ap = ones(1,n+1);
bm = ones(1,m+1); 
am = ones(1,n+1);

%% inicializando parametros a serem adaptados
theta1 = zeros(n-1,1);
theta_n = 0;
theta2 = zeros(n-1,1);
theta_2n = 0;
theta = [theta1; theta_n; theta2; theta_2n];

%% inicializando o regressor
w1 = zeros(n-1,1);
w2 = zeros(n-1,1);
w = [w1; y; w2; r];

%% inicializando os filtros
Af = zeros(n-1,n-1);
Bf = zeros(n-1,1);

%% Definir Planta - editar de acordo com o grau
kp = 2;
%b = 3;
P = planta(ap,bp,kp);

%% Definir Modelo - editar de acordo com o grau
km = 1;
bm = 2;
am = 5;
M = planta(am,bm,km);

%% Definições de Funções
function P = planta(a,b,Kp)
    P =tf(Kp*b, a);
end    