function resultados = comparar_erros_desafio_ls()
% COMPARAR_ERROS_DESAFIO_LS  Compara os erros do desafio usando somente LS.
%
% Coloque este arquivo na pasta do modelo least_square_2ordem_des e execute:
%   resultados_desafio_ls = comparar_erros_desafio_ls;
%
% Este arquivo e independente do script original. Os parametros abaixo
% reproduzem a configuracao fornecida. Nao e necessario rodar o script antigo.
%
% Sao executados tres casos do MESMO modelo:
%   somente_phi1 : gamma1 = 10*I, gamma2 = 0;
%   somente_phi2 : gamma1 = 0,    gamma2 = 10*I;
%   combinado   : gamma1 = 10*I, gamma2 = 10*I.
%
% No modelo do Git, e0 e e sao os erros de PREDICAO dos canais 1 e 2.
% Os erros PARAMETRICOS sao calculados aqui como theta - theta_star.
% Os estados theta1 e theta2 recebem a mesma atualizacao combinada.
% Portanto, theta1 e theta2 nao sao dois estimadores independentes nesse modelo.
%
% As alteracoes de ganhos sao temporarias, via SimulationInput. Este arquivo
% nao salva alteracoes no .mdl/.slx, nao limpa o workspace e nao exporta figuras.
% Mantem o solver e suas tolerancias conforme configurados no modelo local.

modelo = 'least_square_2ordem_des';

%% Parametros do experimento -- mesmos valores do script fornecido
p.tfinal = 500;
p.st = 0.01;                 % Amostragem dos blocos To Workspace.
p.lambda2 = 3;
p.lambda1 = 2;
p.lambda0 = 1;
p.a1 = 2.5;
p.a0 = 1.5;
p.b1 = 2;
p.b0 = 0.8;
p.P = ss(tf([p.b1 p.b0], [1 p.a1 p.a0]));
p.yp0 = [0; 0];
p.x0 = p.yp0;
p.DC = 1;
p.Aq = 2;
p.wq = 0.1*pi;
p.As = 1;
p.ws = pi;
p.theta_star1 = [p.b0; p.b1; p.a0; p.a1];
p.theta_star2 = p.theta_star1;
p.theta0 = zeros(4, 1);
p.kappa = 1;
p.P01 = eye(4);
p.P02 = eye(4);

% Ganhos ativos em cada contribuicao. A comparacao mantem 10*I por canal,
% como no script original, sem dividir o ganho do caso combinado por dois.
ganho1 = 10*eye(4);
ganho2 = 10*eye(4);
zero = zeros(4);

chaves = {'somente_phi1', 'somente_phi2', 'combinado'};
rotulos = {'Somente $\phi_1$', 'Somente $\phi_2$', '$\phi_1+\phi_2$'};
ganhos1 = {ganho1, zero, ganho1};
ganhos2 = {zero, ganho2, ganho2};

%% Simulacoes: salvar os dois canais antes de executar o caso seguinte
load_system(modelo);
resultados = struct('modelo', modelo, 'parametros', p);
nomes = fieldnames(p);

for c = 1:numel(chaves)
    fprintf('LS: executando %s...\n', chaves{c});
    entrada = Simulink.SimulationInput(modelo);
    for k = 1:numel(nomes)
        entrada = entrada.setVariable(nomes{k}, p.(nomes{k}));
    end
    entrada = entrada.setVariable('gamma1', ganhos1{c});
    entrada = entrada.setVariable('gamma2', ganhos2{c});
    entrada = entrada.setModelParameter( ...
        'StartTime', '0', 'StopTime', num2str(p.tfinal, 17), ...
        'LoadInitialState', 'off', 'ReturnWorkspaceOutputs', 'on');
    saida = sim(entrada);

    r = struct;
    r.t = double(lerValor(saida, 't'));
    r.t = r.t(:);
    if numel(r.t) < 2 || any(~isfinite(r.t)) || any(diff(r.t) <= 0)
        error('DesafioLS:TempoInvalido', ...
            'O sinal t deve conter pelo menos duas amostras em ordem crescente.');
    end
    n = numel(r.t);
    r.gamma1 = ganhos1{c};
    r.gamma2 = ganhos2{c};
    r.theta1 = lerSerie(saida, 'theta1', n, 4);
    r.theta2 = lerSerie(saida, 'theta2', n, 4);
    r.phi1 = lerSerie(saida, 'phi1', n, 4);
    r.phi2 = lerSerie(saida, 'phi2', n, 4);

    % Erros de predicao reais, registrados pelo proprio Simulink.
    % Nao usar y-r: esse sinal nao e o erro de predicao da regressao.
    r.e0 = lerSerie(saida, 'e0', n, 1);    % epsilon_1 = theta1'*phi1 - Y1
    r.e  = lerSerie(saida, 'e',  n, 1);    % epsilon_2 = theta2'*phi2 - Y2
    r.erro_predicao = [r.e0 r.e];

    % O vetor identificado e [b0; b1; a0; a1]. Para os graficos usamos
    % theta1, pois ambos os estados recebem a mesma derivada e o mesmo IC.
    r.erro_parametrico = bsxfun(@minus, r.theta1, p.theta_star1.');
    r.erro_parametrico_abs = abs(r.erro_parametrico);
    r.norma_erro_parametrico = sqrt(sum(r.erro_parametrico.^2, 2));
    r.diferenca_entre_estados = max(abs(r.theta1(:) - r.theta2(:)));

    % Verificacao da identidade de erro usada pelo modelo do Git.
    % Sao diagnosticos separados dos sinais reais de predicao acima.
    tilde2 = bsxfun(@minus, r.theta2, p.theta_star2.');
    r.residuo_identidade_predicao = [ ...
        max(abs(r.e0 - sum(r.erro_parametrico .* r.phi1, 2))), ...
        max(abs(r.e  - sum(tilde2 .* r.phi2, 2)))];

    resultados.(chaves{c}) = r;
end

%% Aparencia comum: a cor e o estilo identificam o mesmo caso em cada figura
cores = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.1000 0.1000 0.1000];
estilos = {'--', ':', '-'};
larguras = [1.3 1.5 1.6];
parametros = {'b_0', 'b_1', 'a_0', 'a_1'};

%% Figura 1 -- erro absoluto de cada parametro, com os tres casos sobrepostos
fig = figure('Name', 'Desafio LS - erros parametricos individuais', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [60 70 1200 720]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:4
    ax = nexttile(layout);
    hold(ax, 'on');
    for c = 1:numel(chaves)
        r = resultados.(chaves{c});
        plot(ax, r.t, r.erro_parametrico_abs(:, i), ...
            'Color', cores(c,:), 'LineStyle', estilos{c}, ...
            'LineWidth', larguras(c), 'DisplayName', rotulos{c});
    end
    title(ax, ['$|\widetilde{\theta}_{' num2str(i) '}|$ (' ...
        '$' parametros{i} '$)'], 'Interpreter', 'latex');
    ylabel(ax, '$|\theta_i-\theta_i^*|$', 'Interpreter', 'latex');
    prepararEixo(ax, p.tfinal);
end
title(layout, 'LS: erros parametricos por componente');
resultados.figuras.erros_parametricos = fig;

%% Figura 2 -- norma do erro parametrico para comparar a convergencia global
fig = figure('Name', 'Desafio LS - norma do erro parametrico', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [90 100 1100 450]);
ax = axes('Parent', fig);
hold(ax, 'on');
for c = 1:numel(chaves)
    r = resultados.(chaves{c});
    plot(ax, r.t, r.norma_erro_parametrico, ...
        'Color', cores(c,:), 'LineStyle', estilos{c}, ...
        'LineWidth', larguras(c), 'DisplayName', rotulos{c});
end
ylabel(ax, '$\|\widetilde{\theta}(t)\|_2$', 'Interpreter', 'latex');
title(ax, 'LS: norma do erro parametrico');
prepararEixo(ax, p.tfinal);
resultados.figuras.norma_erro_parametrico = fig;

%% Figuras 3 e 4 -- predicao por canal, com sinal e em valor absoluto
% Mesmo no caso combinado existem DOIS erros de predicao. Cada painel
% compara o MESMO canal nos tres casos. Nao somamos epsilon_1 e epsilon_2.
nomesFig = {'Desafio LS - erros de predicao', ...
            'Desafio LS - erros absolutos de predicao'};
titulosFig = {'LS: erros de predicao por canal', ...
              'LS: valores absolutos dos erros de predicao'};
titulosCanal = {'Canal 1: sinal e0 do modelo', 'Canal 2: sinal e do modelo'};

for absoluto = 0:1
    fig = figure('Name', nomesFig{absoluto + 1}, ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [120 80 1100 680]);
    layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for canal = 1:2
        ax = nexttile(layout);
        hold(ax, 'on');
        for c = 1:numel(chaves)
            r = resultados.(chaves{c});
            dados = r.erro_predicao(:, canal);
            if absoluto
                dados = abs(dados);
            end
            plot(ax, r.t, dados, 'Color', cores(c,:), ...
                'LineStyle', estilos{c}, 'LineWidth', larguras(c), ...
                'DisplayName', rotulos{c});
        end
        title(ax, titulosCanal{canal});
        if absoluto
            ylabel(ax, ['$|\epsilon_' num2str(canal) '|$'], 'Interpreter', 'latex');
        else
            ylabel(ax, ['$\epsilon_' num2str(canal) '$'], 'Interpreter', 'latex');
        end
        prepararEixo(ax, p.tfinal);
    end
    title(layout, titulosFig{absoluto + 1});
    if absoluto
        resultados.figuras.predicao_absoluta = fig;
    else
        resultados.figuras.predicao = fig;
    end
end

%% Resumo numerico e acesso aos dados
fprintf('\nNorma do erro parametrico no ultimo instante:\n');
for c = 1:numel(chaves)
    r = resultados.(chaves{c});
    fprintf('  %-14s : %.6g\n', chaves{c}, r.norma_erro_parametrico(end));
end
fprintf('\ne0 = erro de predicao do canal 1; e = erro de predicao do canal 2.\n');
fprintf('Erros parametricos assinados: resultados.<caso>.erro_parametrico\n');
fprintf('Erros de predicao: resultados.<caso>.e0 e resultados.<caso>.e\n');
end

function prepararEixo(ax, tfinal)
grid(ax, 'on');
xlabel(ax, 'Tempo (s)');
xlim(ax, [0 tfinal]);
set(ax, 'FontSize', 11);
legend(ax, 'show', 'Interpreter', 'latex', 'Location', 'best');
end

function valor = lerValor(saida, nome)
% Os nomes abaixo sao os To Workspace existentes no modelo do Git.
if ~ismember(nome, who(saida))
    error('DesafioLS:SinalAusente', ...
        ['O modelo nao retornou o sinal "%s". Verifique o bloco To Workspace ' ...
         'com esse VariableName no modelo least_square_2ordem_des.'], nome);
end
valor = saida.get(nome);
if isa(valor, 'Simulink.SimulationData.Signal')
    valor = valor.Values;
end
if isa(valor, 'timeseries')
    valor = valor.Data;
elseif isstruct(valor) && isfield(valor, 'signals')
    valor = valor.signals.values;
end
if ~isnumeric(valor) || ~isreal(valor)
    error('DesafioLS:FormatoInvalido', ...
        'O sinal %s deve conter um array numerico real.', nome);
end
end

function dados = lerSerie(saida, nome, nAmostras, nCanais)
% Aceita arrays Nxd, dxN e logs de vetores em formato dx1xN.
dados = double(squeeze(lerValor(saida, nome)));
if nCanais == 1 && isvector(dados)
    dados = dados(:);
elseif isequal(size(dados), [nCanais nAmostras])
    dados = dados.';
end
if ~isequal(size(dados), [nAmostras nCanais])
    error('DesafioLS:DimensaoInvalida', ...
        'Sinal %s: esperado %dx%d, recebido %s.', ...
        nome, nAmostras, nCanais, mat2str(size(dados)));
end
if any(~isfinite(dados(:)))
    error('DesafioLS:SinalNaoFinito', ...
        'O sinal %s possui NaN ou Inf. Verifique a simulacao.', nome);
end
end
