function SolucaoOtimizador = CriaSolver_NMPCESNx(Hp,Hc,q,r,qu,ny,nu,nx,ModeloPreditor)
% InicalizaSolver_NMPC ->Cria o otimizador para NMPC.
%
% Entradas:
%   Hp   - Horizonte de predição.
%   Hc   - Horizonte de controle.
%   q    - Vetor de ponderação para o erro de rastreamento da saída.
%   r    - Vetor de ponderação para os incrementos da entrada de controle.
%   qu   - Vetor de ponderação para o desvio da entrada de controle em relação a um alvo.
%   ny   - Número de saídas controladas.
%   nu   - Número de entradas de controle.
%   nx   - Número de estados.
%   ModeloPreditor - Modelo da rede ESN identificada pelo Python
%
% Saída:
%   SolucaoOtimizador - Objeto nlpsol do CasADi para o problema NMPC.
import casadi.*
%% ========================================================================
%  1. FORMULAÇÃO DO PROBLEMA NMPC
%  ========================================================================

% ------------------------------------------------------------------------
%  1.1. Define Parâmetros do NMPC, Matrizes de Ponderação e Variáveis Simbólicas
% ------------------------------------------------------------------------
Qy = diag(q);                % Matriz de ponderação para o erro de rastreamento da saída
R = diag(repmat(r,1,Hc));    % Matriz de ponderação para os incrementos da entrada de controle (sobre Hc)
Qu = diag(qu);               % Matriz de ponderação para o desvio da entrada em relação a um alvo

% Variáveis simbólicas para o problema de otimização
X =     MX.sym('X',nx,Hp+1);                % Trajetória de estados predita sobre Hp
du =    MX.sym('du',Hc*nu,1);               % Variável de decisão: sequência de incrementos da entrada de controle sobre Hc
Du =    [du;zeros(nu*(Hp-Hc),1)];           % Sequência completa de incrementos da entrada (preenchida para Hp > Hc)
ysp =   MX.sym('ysp',ny,1);                 % Variável de decisão: setpoint (se otimizado, ou um parâmetro se fixo)

% Define variáveis controladas como uma função dos estados
VarControladas = MX.sym('VarControladas',nx);          
h_ESN=Function('h_ESN',{VarControladas},  {[VarControladas(1);VarControladas(2)]});  % cria o vetor das medições

% ------------------------------------------------------------------------
%  1.2. Define Parâmetros Simbólicos de Tempo de Execução para o Solver
% ------------------------------------------------------------------------
nx_ESN = length(ModeloPreditor.data.a0);     % Define o tamanho da rede ESN     
P =      MX.sym('P',nx+nu+ny+nu+nu+nx_ESN);  % P é um vetor de parâmetros externos passados ao solver a cada chamada.
uk_1 =   P(ny+nu:ny+nu+1);                   % Ação de controle anterior
er=      P(ny+nu+nu:ny+nu+nu+1);             % Erro anterior da medição entre a planta e o modelo preditor
uRTO =   P(ny+nu+nu+nu+nu:ny+nu+nu+nu+nu+1); % Target desejado
ESNdataa0 =   P(ny+nu+nu+nu+nu+nu:end);      % Estados do reservatório

% ------------------------------------------------------------------------
%  1.3. Formula Função Objetivo (fob) e Restrições (g)
% ------------------------------------------------------------------------
fob = 0; % Inicializa função objetivo
g = [];  % Inicializa vetor de restrições
g=[g;X(:,1)-P(1:nx)];  % Restrição inicial: Estado predito em k=0 é o estado medido atual
[u_max,u_min,y_max,y_min] = normalizar_BCS_ESN; % Obtém parâmetros de normalização da ESN 

% Loop sobre o horizonte de predição Hp
for k=1:Hp                                  
        % Dentro do horizonte de controle (Hc), usa incrementos 'du'. Após u permanece constante.
    if k<Hc
        duk_aux = Du((k-1)*nu+1:k*nu);
        uk_1 = uk_1 + duk_aux;         
    end
    ym = h_ESN(X(:,k+1));     % Saída predita no próximo passo (k+1)
    fob=fob + (ym-ysp+er)'*Qy*(ym-ysp+er)+du'*R*du+(uk_1-uRTO)'*Qu*(uk_1-uRTO);            % define a função objetivo proposta
    %Atualiza o modelo preditor com os parametros de entrada e regressores 
    ukk_aux = feature_scaling(vertcat(uk_1,P(8:9),P(1:3)),u_max(1:7),u_min(1:7));
    x_ESN = ModeloPreditor.data.Wrr*ESNdataa0 + ModeloPreditor.data.Wir*ukk_aux + ModeloPreditor.data.Wbr;  %usar o modeloPreditor(ESN) para fazer a predição
    next_state = (1-ModeloPreditor.data.gama)*ESNdataa0 + ModeloPreditor.data.gama*tanh(x_ESN);         % ||
    a_wbias = [1.0;next_state];                                                                         % ||
    yn = ModeloPreditor.data.Wro*a_wbias;  
    y_esn_pred = feature_descaling(yn,y_max(1:3),y_min(1:3));
    
    g=[g;X(:,k+1)-y_esn_pred]; % Restrição de igualdade: X(k+1) - X_next_pred = 0
end
% ------------------------------------------------------------------------
%  1.4. Define Matrizes Auxiliares (Mtil, Itil) para Restrições de Entrada
% ------------------------------------------------------------------------
% Estas matrizes são usadas para formular restrições nas entradas de controle reais uk
Mtil=[];                         
Itil=[];
auxM=zeros(nu,Hc*nu);
for in=1:Hc
    auxM=[eye(nu) auxM(:,1:(Hc-1)*nu)];
    Mtil=[Mtil;auxM];
    Itil=[Itil;eye(nu)];
end
% Conclui a inclusão das restrições nas entradas para os limites: u_min <= Mtil*du + Itil*uk <= u_max
g = [g;Mtil*du+Itil*P(nx+1:nx+nu)]; 

%% ========================================================================
%  2. CONFIGURAÇÃO DO SOLVER NLP (Programação Não Linear)
%  ========================================================================

% ------------------------------------------------------------------------
%  2.1. Define Variáveis de Otimização e Estrutura do Problema NLP
% ------------------------------------------------------------------------
% Variáveis a serem otimizadas pelo solver
opt_variable=[X(:);du;ysp];  

% Estrutura do problema NLP para o CasADi
% 'f': função objetivo
% 'x': variáveis de decisão (opt_variables)
% 'g': funções de restrição
% 'p': parâmetros de tempo de execução (P)
nlp = struct('f',fob,'x',opt_variable,'g', g, 'p', P); %define a estrutura para problema de otimização não linear (NLP, Nonlinear Programming), sendo: 1-fob, definida acima; 2-opt_variable: variáveis de decisão; 3-g, as restrições do processo e 4-P, parâmetros de entrada para Solver   

% ------------------------------------------------------------------------
%  2.2. Define Opções do Solver (IPOPT)
% ------------------------------------------------------------------------
options=struct;
options.print_time=0;                         % Habilita tempo total de execução do solver deve ser impresso ou não.
options.ipopt.print_level=1;                  % Nível de detalhamento das mensagens de saída do IPOPT. Valores mais baixos resultam em menos mensagens (0 significa sem mensagens).
options.ipopt.max_iter=1000;                  % Especifica o número máximo de iterações que o solver deve executar antes de parar.
options.ipopt.acceptable_tol=1e-8;            % Define a tolerância de convergência do solver. Um valor menor indica uma solução mais precisa.
options.ipopt.acceptable_obj_change_tol=1e-8; % Define uma tolerância aceitável para uma solução "boa o suficiente", útil para problemas onde a solução perfeita pode ser muito difícil de alcançar.
options.ipopt.max_wall_time=30;               % Tempo (em segundos) máximo para o solver encontrar solução
% ------------------------------------------------------------------------
%  2.3. Inicializa o Solver NLP
% ------------------------------------------------------------------------
SolucaoOtimizador = nlpsol('SolucaoOtimizador','ipopt', nlp,options); 
