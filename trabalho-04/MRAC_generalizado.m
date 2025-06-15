%% graus da planta
n = 5; %editar aqui - grau do denominador
m = 3; %editar aqui - grau do numerador
n_star = n-m; %editar aqui - g1rau relativo


%% inicializando variaveis
kp = 0;
km = 0;
bp = zeros(1,m+1); 
ap = zeros(1,n+1);
bm = zeros(1,m+1); 
am = zeros(1,n+1);
theta = zeros(2*n);
w = zeros(2*n);
w1 = zeros(n-1);
w2 = zeros(n-1);
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




%%
function P = planta(a,b,Kp)
    P =tf(Kp*b, a);
end    