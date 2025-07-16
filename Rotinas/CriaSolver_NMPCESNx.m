function SolucaoOtimizador = CriaSolver_NMPCESNx(Hp,Hc,q,r,qu,ny,nu,nx,ModeloPreditor)
% InicalizaSolver_NMPC ->Cria o otimizador para NMPC.
%
% Entradas:
%   Hp   - Horizonte de predi��o.
%   Hc   - Horizonte de controle.
%   q    - Vetor de pondera��o para o erro de rastreamento da sa�da.
%   r    - Vetor de pondera��o para os incrementos da entrada de controle.
%   qu   - Vetor de pondera��o para o desvio da entrada de controle em rela��o a um alvo.
%   ny   - N�mero de sa�das controladas.
%   nu   - N�mero de entradas de controle.
%   nx   - N�mero de estados.
%   ModeloPreditor - Modelo da rede ESN identificada pelo Python
%
% Sa�da:
%   SolucaoOtimizador - Objeto nlpsol do CasADi para o problema NMPC.
import casadi.*
%% ========================================================================
%  1. FORMULA��O DO PROBLEMA NMPC
%  ========================================================================

% ------------------------------------------------------------------------
%  1.1. Define Par�metros do NMPC, Matrizes de Pondera��o e Vari�veis Simb�licas
% ------------------------------------------------------------------------
Qy = diag(q);                % Matriz de pondera��o para o erro de rastreamento da sa�da
R = diag(repmat(r,1,Hc));    % Matriz de pondera��o para os incrementos da entrada de controle (sobre Hc)
Qu = diag(qu);               % Matriz de pondera��o para o desvio da entrada em rela��o a um alvo

% Vari�veis simb�licas para o problema de otimiza��o
X =     MX.sym('X',nx,Hp+1);                % Trajet�ria de estados predita sobre Hp
du =    MX.sym('du',Hc*nu,1);               % Vari�vel de decis�o: sequ�ncia de incrementos da entrada de controle sobre Hc
Du =    [du;zeros(nu*(Hp-Hc),1)];           % Sequ�ncia completa de incrementos da entrada (preenchida para Hp > Hc)
ysp =   MX.sym('ysp',ny,1);                 % Vari�vel de decis�o: setpoint (se otimizado, ou um par�metro se fixo)

% Define vari�veis controladas como uma fun��o dos estados
VarControladas = MX.sym('VarControladas',nx);          
h_ESN=Function('h_ESN',{VarControladas},  {[VarControladas(1);VarControladas(2)]});  % cria o vetor das medi��es

% ------------------------------------------------------------------------
%  1.2. Define Par�metros Simb�licos de Tempo de Execu��o para o Solver
% ------------------------------------------------------------------------
nx_ESN = length(ModeloPreditor.data.a0);     % Define o tamanho da rede ESN     
P =      MX.sym('P',nx+nu+ny+nu+nu+nx_ESN);  % P � um vetor de par�metros externos passados ao solver a cada chamada.
uk_1 =   P(ny+nu:ny+nu+1);                   % A��o de controle anterior
er=      P(ny+nu+nu:ny+nu+nu+1);             % Erro anterior da medi��o entre a planta e o modelo preditor
uRTO =   P(ny+nu+nu+nu+nu:ny+nu+nu+nu+nu+1); % Target desejado
ESNdataa0 =   P(ny+nu+nu+nu+nu+nu:end);      % Estados do reservat�rio

% ------------------------------------------------------------------------
%  1.3. Formula Fun��o Objetivo (fob) e Restri��es (g)
% ------------------------------------------------------------------------
fob = 0; % Inicializa fun��o objetivo
g = [];  % Inicializa vetor de restri��es
g=[g;X(:,1)-P(1:nx)];  % Restri��o inicial: Estado predito em k=0 � o estado medido atual
[u_max,u_min,y_max,y_min] = normalizar_BCS_ESN; % Obt�m par�metros de normaliza��o da ESN 

% Loop sobre o horizonte de predi��o Hp
for k=1:Hp                                  
        % Dentro do horizonte de controle (Hc), usa incrementos 'du'. Ap�s u permanece constante.
    if k<Hc
        duk_aux = Du((k-1)*nu+1:k*nu);
        uk_1 = uk_1 + duk_aux;         
    end
    ym = h_ESN(X(:,k+1));     % Sa�da predita no pr�ximo passo (k+1)
    fob=fob + (ym-ysp+er)'*Qy*(ym-ysp+er)+du'*R*du+(uk_1-uRTO)'*Qu*(uk_1-uRTO);            % define a fun��o objetivo proposta
    %Atualiza o modelo preditor com os parametros de entrada e regressores 
    ukk_aux = feature_scaling(vertcat(uk_1,P(8:9),P(1:3)),u_max(1:7),u_min(1:7));
    x_ESN = ModeloPreditor.data.Wrr*ESNdataa0 + ModeloPreditor.data.Wir*ukk_aux + ModeloPreditor.data.Wbr;  %usar o modeloPreditor(ESN) para fazer a predi��o
    next_state = (1-ModeloPreditor.data.gama)*ESNdataa0 + ModeloPreditor.data.gama*tanh(x_ESN);         % ||
    a_wbias = [1.0;next_state];                                                                         % ||
    yn = ModeloPreditor.data.Wro*a_wbias;  
    y_esn_pred = feature_descaling(yn,y_max(1:3),y_min(1:3));
    
    g=[g;X(:,k+1)-y_esn_pred]; % Restri��o de igualdade: X(k+1) - X_next_pred = 0
end
% ------------------------------------------------------------------------
%  1.4. Define Matrizes Auxiliares (Mtil, Itil) para Restri��es de Entrada
% ------------------------------------------------------------------------
% Estas matrizes s�o usadas para formular restri��es nas entradas de controle reais uk
Mtil=[];                         
Itil=[];
auxM=zeros(nu,Hc*nu);
for in=1:Hc
    auxM=[eye(nu) auxM(:,1:(Hc-1)*nu)];
    Mtil=[Mtil;auxM];
    Itil=[Itil;eye(nu)];
end
% Conclui a inclus�o das restri��es nas entradas para os limites: u_min <= Mtil*du + Itil*uk <= u_max
g = [g;Mtil*du+Itil*P(nx+1:nx+nu)]; 

%% ========================================================================
%  2. CONFIGURA��O DO SOLVER NLP (Programa��o N�o Linear)
%  ========================================================================

% ------------------------------------------------------------------------
%  2.1. Define Vari�veis de Otimiza��o e Estrutura do Problema NLP
% ------------------------------------------------------------------------
% Vari�veis a serem otimizadas pelo solver
opt_variable=[X(:);du;ysp];  

% Estrutura do problema NLP para o CasADi
% 'f': fun��o objetivo
% 'x': vari�veis de decis�o (opt_variables)
% 'g': fun��es de restri��o
% 'p': par�metros de tempo de execu��o (P)
nlp = struct('f',fob,'x',opt_variable,'g', g, 'p', P); %define a estrutura para problema de otimiza��o n�o linear (NLP, Nonlinear Programming), sendo: 1-fob, definida acima; 2-opt_variable: vari�veis de decis�o; 3-g, as restri��es do processo e 4-P, par�metros de entrada para Solver   

% ------------------------------------------------------------------------
%  2.2. Define Op��es do Solver (IPOPT)
% ------------------------------------------------------------------------
options=struct;
options.print_time=0;                         % Habilita tempo total de execu��o do solver deve ser impresso ou n�o.
options.ipopt.print_level=1;                  % N�vel de detalhamento das mensagens de sa�da do IPOPT. Valores mais baixos resultam em menos mensagens (0 significa sem mensagens).
options.ipopt.max_iter=1000;                  % Especifica o n�mero m�ximo de itera��es que o solver deve executar antes de parar.
options.ipopt.acceptable_tol=1e-8;            % Define a toler�ncia de converg�ncia do solver. Um valor menor indica uma solu��o mais precisa.
options.ipopt.acceptable_obj_change_tol=1e-8; % Define uma toler�ncia aceit�vel para uma solu��o "boa o suficiente", �til para problemas onde a solu��o perfeita pode ser muito dif�cil de alcan�ar.
options.ipopt.max_wall_time=30;               % Tempo (em segundos) m�ximo para o solver encontrar solu��o
% ------------------------------------------------------------------------
%  2.3. Inicializa o Solver NLP
% ------------------------------------------------------------------------
SolucaoOtimizador = nlpsol('SolucaoOtimizador','ipopt', nlp,options); 
