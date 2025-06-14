function [theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0)
    % ===============================
    % Função para cálculo dos vetores θ* (controle 2DOF)
    % Entrada:
    %  - P : função de transferência da planta
    %  - M : função de transferência do modelo de referência
    %  - A0: polinômio característico do observador (como vetor)
    % Saída:
    %  - theta1 : vetor do filtro na entrada
    %  - theta_n: ganho direto na saída
    %  - theta2 : vetor do filtro na saída
    %  - theta_2n: ganho no feedforward (referência)
    % ===============================
    
    %% Cálculo simbólico
    syms s
    
    % ===============================
    % Obtenção dos dados da planta e modelo
    % ===============================
    [N_p, D_p] = tfdata(P, 'v');
    [N_m, D_m] = tfdata(M, 'v');

    % O kp está no termo de maior grau (s^1)
    kp = N_p(find(N_p ~= 0, 1, 'first'));
    km = N_m(find(N_m ~= 0, 1, 'first'));
    Np = poly2sym(N_p, s)/kp;
    Dp = poly2sym(D_p, s);

    Nm = poly2sym(N_m, s)/km;
    Dm = poly2sym(D_m, s);

    A0_sym = poly2sym(A0, s);
    
    % ===============================
    % Graus dos polinômios
    % ===============================
    n = length(D_p) - 1;
    m = length(N_p(find(N_p~=0, 1, 'first'):end)) - 1;
    n_star = n - m;

    %% ===============================
    % Cálculo de Λ(s)
    % ===============================
    Lambda = Nm * A0_sym;

    %% ===============================
    % Resolver equação Diofantina:
    % H(s)*D_p + G(s) = Dm * A0
    % ===============================
    deg_H = n_star - 1;
    deg_G = n - 1;

    H = sym('h', [1 deg_H]);
    G = sym('g', [1 deg_G + 1]);

    H_poly = poly2sym([1 H], s);
    G_poly = poly2sym(G, s);

    % Montagem da equação Diofantina
    Eq = expand(H_poly * Dp - kp*G_poly - Dm * A0_sym);

    % Coletar os coeficientes em potências de s
    coeffs_Eq = fliplr(coeffs(Eq, s));

    % Sistema de equações
    eqns = coeffs_Eq == zeros(1, length(coeffs_Eq));

    % Resolver equações
    vars = [H G];
    sol = solve(eqns, vars);

    % Obter nomes dos campos da solução
    fields = fieldnames(sol);

    %% Processar H
    h_fields = fields(startsWith(fields, 'h'));

    if isempty(h_fields)
        H_sol = [];  % Caso não exista H (grau zero)
    else
        h_fields = sort(h_fields); % Ordena lexicograficamente
        H_cells = cellfun(@(c) sol.(c), h_fields, 'UniformOutput', false);
        H_sol = double([H_cells{:}]); % Vetor linha com os coeficientes
    end

    %% Processar G
    g_fields = fields(startsWith(fields, 'g'));

    if isempty(g_fields)
        G_sol = [];  % Caso não exista G
    else
        g_fields = sort(g_fields); % Ordena lexicograficamente
        G_cells = cellfun(@(c) sol.(c), g_fields, 'UniformOutput', false);
        G_sol = double([G_cells{:}]); % Vetor linha com os coeficientes
    end

    
    % Montar os polinômios H(s) e G(s)
    H_final = poly2sym([1 H_sol], s);
    G_final = poly2sym(G_sol, s);

    %% ===============================
    % Cálculo de F(s)
    % F = Λ - N_p * H
    % ===============================
    F = expand(Lambda - Np * H_final);

    %% ===============================
    % Extração dos vetores θ*
    % ===============================

    % θ_2n = km / kp
    % Obter ganhos DC da planta e do modelo
    theta_2n = km / kp;

    % θ2 -> Coeficientes de G(s) normalizados por Λ
    % θn -> ganho direto (coeficiente de y)
    % Divisão polinomial simbólica: G = q*Lambda + r
    [q, r] = quorem(G_final, Lambda, s);

    % θ_n é o quociente (termo sem denominador)
    theta_n = sym2poly(q);

    % θ_2 é o resto (termo que divide Lambda)
    theta2 = sym2poly(r);

    % θ1 -> Coeficientes de F(s) normalizados por Λ
    theta1 = sym2poly(F);

    %% ===============================
    % Saída formatada
    % ===============================
    fprintf('\n=== Resultados ===\n');
    fprintf('theta_1 = %s\n', num2str(theta1));
    fprintf('theta_n = %s\n', num2str(theta_n));
    fprintf('theta_2 = %s\n', num2str(theta2));
    fprintf('theta_2n = %s\n', num2str(theta_2n));
end
