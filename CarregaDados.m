%% Script Principal: Implementação do NMPC orientado a dados no LEA-BCS
% Este script implementa o NMPC+ESNx no ESP-LEA em ambiente de simulação
% Autor: Odilon Abreu
% Data: 1/10/2024

clc         % Limpa a janela de comando.
clear       % Limpa todas as variáveis do workspace.
close all   % Fecha todas as figuras abertas.

% Aponta para diretórios de uso geral
addpath('ModeloESNx');       % Dados do BCS
addpath('Rotinas');     % Rotinas Uteis do BCS

%% 1. Parâmetros Gerais de Simulação
% Define parâmetros fundamentais para a simulação do sistema e do NMPC.
Ts      = 1;                  % Período de amostragem do processo (s).
Ts_NMPC = 30;                 % Passo de amostragem do NMPC.
Par     = dados_LEA_Par(Ts);  % Carrega parâmetros do LEA-BCS para o modelo não linear.


%% 2. Configuração do Cenário e Parâmetros do Processo
% Define as condições operacionais e carrega os parâmetros do modelo LEA-BCS.
NivelRuidoHead = 2;     % Nível de ruído aditivo para a medição de Head da bomba (0 a 10).
NivelRuidoLa   = 3;     % Nível de ruído aditivo para a medição do Nível do anular (0 a 10).
NivelRuido     = 2;     % Nível de ruído geral (pode ser ajustado para outras variáveis).

%% 3. Definição dos Cenários de Operação para o Processo
    Time_Stop    = 3.5 * 3600 * Ts;     % Tempo total de simulação (s).
    La_x0        = 7.42;                % (m) Nível inicial do anular.
    pwh_x0       = 8.91 * 1e05;         % (Pa) Pressão inicial na choke.
    qm_x0        = 1.87/3600;           % (m³/s) Vazão inicial na coluna de produção.
    head_x0      = 111;                 % (m) Head inicial da bomba.
    ukk          = [50;25];             % (Hz e %) Condições iniciais das entradas (frequência e abertura da choke).
    La_SP        = 6;                 % Setpoint de nível do anular.

   % Perfil de Pressão (reservatório e manifold)  
    Par.Resv = ones(Time_Stop,1)* 1.56; % Perfil de Pressão do Reservatório
    Par.Pman = ones(Time_Stop,1)* 3.0;  % Perfil de Pressão manifold.                       
    Par.Tintake = ones(Time_Stop,1)* Par.Tintake; % Perfil de Temperatura na entrada da Bomba.
    
    % Perfil de controle manual (frequência e choke) para estabilização inicial.
    tempo_manual = 450; % Instante inicial para manter a operação estacionária (s).
    Zc_manual    = ones(1, tempo_manual) * ukk(2);
    f_manual     = ones(1, tempo_manual) * ukk(1);
    tempo_target = 11000;                % Instante para acionar o target de Zc.

%% 4. Parâmetros do modelo não linear (E.D.O)
% Carrega o modelo não linear que representa o PROCESSO (LEA-BCS) e define
% os estados iniciais detalhados do sistema.
dinamica_BCS; % Script que define a dinâmica não linear do LEA-BCS.
Pin        = 7.3e+05;       % (Pa) Pressão de intake (entrada da bomba).
dp_bomba   = 6.33e+05;      % (Pa) Delta de pressão da bomba.
P_hp_bomba = 1.07;          % (HP) Potência requerida pela bomba.
vazao_choke= qm_x0;         % (m³/s) Vazão através da choke.
qr_res     = qm_x0;         % (m³/s) Vazão no reservatório.

% Vetor de estados iniciais do processo para simulação.
xmk = [La_x0; pwh_x0; qm_x0; Pin; head_x0; dp_bomba; P_hp_bomba; vazao_choke; qr_res; Par.visc; Par.IP];

%% 5. Parâmetros de Sintonia do NMPC
% Define os horizontes e pesos para a função objetivo do NMPC.
Hp      = 10;               % Horizonte de predição.
Hc      = 3;                % Horizonte de controle.
q       = [1e2 10];         % Matriz de peso para as saídas (erros). [La, Head].
r       = [0.5,1];          % Matriz de supressão para o incremento das entradas (freq, Zc).
qu      = [0 1];            % Matriz de peso para o desvio das entradas em relação ao target (otimização econômica).
dumax   = [2  1]';          % Variação máxima permitida para as ações de controle (dFreq, dZc).
udes =    [0 46]';          % Definição do Target

du0     = zeros(Hc * size(dumax,1), 1); % Estimativa inicial para os incrementos de controle.
YIni=  [La_x0;head_x0;qm_x0*3600];
X0_NMPC = repmat(YIni,1,Hp+1); % Estimativa inicial para os estados preditos.
% Estimativa inicial para o solver do NMPC (variáveis de decisão empilhadas).
X0Du0YsP = [X0_NMPC(:);du0;[La_SP*1.05;head_x0]]; % [Estados_preditos; Incrementos_de_controle; Setpoint_de_saída].

%% 6. Carregamento ESNx
NomeESN='weightsESNx_ESPLEA.mat';
ESN = load(NomeESN);
ModeloPreditor = ESN;
disp(['Dados carregados para ' NomeESN]);

%% 8. Fim da Configuração
disp('Configuração concluída. Pronto para executar o "NMPC_ESNx_ESPLEA.slx"');



