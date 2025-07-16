classdef casadi_blockNMPC_ESNx < matlab.System & matlab.system.mixin.Propagates
    % casadi_blockNMPC_ESNx - Bloco MATLAB System para implementação do NMPC com CasADi.
    % Este bloco gerencia a lógica de execução do NMPC no Simulink,
    % incluindo a inicialização do solver, manipulação de dados de entrada
    % e saída, aplicação de restrições e gestão de cenários de operação.
    properties
    end
    
    properties (DiscreteState)
    end

    properties (Access = private)
        % Propriedades privadas para armazenamento interno do estado do bloco.
        casadi_solver           % Solver CasADi para NMPC com objetivo econômico.
        x0                      % Estimativa inicial para as variáveis de decisão do solver.
        uk                      % Ação de controle atual (frequência e abertura da choke)
        contador                % Contador para o PassoMPC (intervalo de execução do NMPC)
        Tsolver                 % Tempo de cálculo da última execução do solver
        Envelope                % Estrutura carregada do arquivo .mat com restrições de envelope
        fob                     % Valor da função objetivo da última otimização
        ysp                     % Setpoint(s) [ysp1; ysp2] (pode ser otimizado pelo NMPC)
        ModeloPreditor          % Estrutura ou objeto contendo o modelo preditivo (ex: ESN)
        f_manual                % Vetor de frequência manual para controle inicial
        Zc_manual               % Vetor de abertura de choke manual para controle inicial
        tempo_manual            % Tempo durante o qual o controle manual é aplicado
        tempo_target            % Tempo para alguma lógica de alvo (não usado explicitamente no stepImpl fornecido)
        duk_anterior
    end

    methods (Access = protected)
        %% 1. Métodos de Configuração de Portas de Entrada/Saída (para Simulink)
        %   Estes métodos informam ao Simulink sobre as características das
        %   portas do bloco (número, tipo de dado, tamanho, complexidade).

        function num = getNumInputsImpl(~)
             num = 11; % 4 (DadosProcesso assume-se dividido) + 6 (escalares) + 1 (PassoMPC) -> Ver stepImpl para as 11 entradas
                      % Ou: 1 (DadosProcesso) + 1 (t) + 1 (SP) + 1 (uRTO) + 1 (dumax) + 1 (Hp) + 1 (Hc) + 1 (q) + 1 (r) + 1 (qu) + 1 (PassoMPC) = 11
        end

        function num = getNumOutputsImpl(~)
        % Define o número de portas de saída do bloco.
            num = 1; % Uma única porta de saída, mas é um vetor grande.
        end

        function dt1 = getOutputDataTypeImpl(~)
        % Define o tipo de dados da porta de saída.
            dt1 = 'double';
        end

        function dt1 = getInputDataTypeImpl(~)
        % Define o tipo de dados das portas de entrada. (Assumido 'double' para todas)
            dt1 = 'double';
        end

        function sz1 = getOutputSizeImpl(~)
        % Define o tamanho da porta de saída.
        % Saída é um vetor concatenado: [uk(2); ysp(2); ymax(2); ymin(2); LimitesMax(4, inc H_max); LimitesMin(4, inc H_min); ManipHigh(2); ManipLow(2); fob(1); Tsolver(1); Feasible(1); PredHp(3); PredHead(1)]
            sz1 = [2+2+2+2+4+4+2+2+1+1+1+3+1]; % Total 27 elementos
        end

        function sz1 = getInputSizeImpl(~)
        % Define o tamanho de cada porta de entrada. 
            sz1 = [1,8,1,2,2,1,1,1,1,1,1,1]; % Esta definição de tamanho pode ser para um cenário específico.
        end

        function cp1 = isInputComplexImpl(~)
        % Define se as entradas são números complexos.
            cp1 = false;
        end

        function cp1 = isOutputComplexImpl(~)
        % Define se as saídas são números complexos.
            cp1 = false;
        end

        function fz1 = isInputFixedSizeImpl(~)
        % Define se o tamanho das entradas é fixo.
            fz1 = true;
        end

        function fz1 = isOutputFixedSizeImpl(~)
        % Define se o tamanho das saídas é fixo.
            fz1 = true;
        end
        
        %% 2. Método de Inicialização do Bloco (setupImpl)
        function setupImpl(obj,~)
            % SETUPIMPL - Executado uma única vez no início da simulação.
            obj.casadi_solver = [];           % Inicializa o solver com target econômico.

            % Carrega variáveis do workspace base do MATLAB (definidas no script principal).
            obj.x0 =  evalin('base','X0Du0YsP');      % Palpite inicial para o solver NMPC
            obj.uk = evalin('base', 'ukk');           % Última ação de controle conhecida/inicial
            obj.ModeloPreditor = evalin('base', 'ModeloPreditor'); % Carrega o modelo preditor (ex: ESN)
            obj.f_manual = evalin('base', 'f_manual');             % Carrega sequência de frequência manual
            obj.Zc_manual = evalin('base', 'Zc_manual');           % Carrega sequência de choke manual
            obj.tempo_manual = evalin('base', 'tempo_manual');     % Carrega tempo de controle manual
            obj.tempo_target = evalin('base', 'tempo_target');     % Carrega tempo de alvo (uso não claro aqui)
            
            obj.contador = 0;                                      % Zera o contador para o PassoMPC
            obj.Tsolver = 0;                                       % Zera o tempo de cômputo do solver

            NomeEnvelope='restricao_Envelope.mat'; % Nome do arquivo de restrições
            obj.Envelope = load(NomeEnvelope);        % Carrega o arquivo de envelope
            
            obj.fob = 0;                  % Inicializa a função objetivo.
            obj.ysp = obj.x0(end-1:end);  % Setpoint inicial (últimos 2 elementos de X0Du0YsP). 

        end
        
        %% 3. Método Principal de Execução (stepImpl)
        function u = stepImpl(obj,DadosProcesso,t,SP,uRTO,dumax,Hp,Hc,q,r,qu,PassoMPC)
            % STEPIMPL - Executado a cada passo de tempo da simulação.
            % Calcula a ação de controle do NMPC e outras saídas.
            
%             disp(strcat("Simulacao MPC_ESN_LEABCS em ",num2str(t)," seg ")); % Mensagem de log
            
            % Dimensões das variáveis do sistema.
            ny = size(q', 2); % Número de saídas (yk).
            nu = size(r', 2); % Número de entradas manipuladas (uk).
            nx = size(DadosProcesso(1:3),1);% Número de estados (La, Pwh, qm).
            
            % Obtém parâmetros de normalização da ESN
            [u_max,u_min,~,~] = normalizar_BCS_ESN;

            % Bloco de inicialização do solver (executado apenas em t=0)
            if t==0
                obj.casadi_solver = CriaSolver_NMPCESNx(Hp,Hc,q',r',qu,ny,nu,nx,obj.ModeloPreditor); % cria o solver (otimizador) uma vez
                % Prepara dados de entrada para o "aquecimento" da ESN (inclui estados)
                ukk = feature_scaling(vertcat(obj.uk,DadosProcesso(4:5),DadosProcesso(1:3)),u_max(1:7),u_min(1:7));
                % "Aquece" a ESN: atualiza o estado do reservatório da ESN várias vezes para estabilização inicial
                obj.ModeloPreditor.data.a0 = esquenta_ESNx(obj.ModeloPreditor.data,ukk,1000); % Atualiza várias vezes o estado do reservátório para esquentar ESN
            end
            
            %% 3.1. Definição dos Limites de Variáveis Manipuladas (Controle)
            umax = [60; 100]; % Limites superiores para frequência (Hz) e choke (%).
            umin = [40; 1];   % Limites inferiores para frequência (Hz) e choke (%).
            % Estes são os limites 'lbg' e 'ubg' do solver.
            ManipuladasLowLimit  = [repmat(zeros(nx,1), Hp+1, 1); repmat(umin, Hc, 1)]; % Limites inferiores de g para X e u 
            ManipuladasHighLimit = [repmat(zeros(nx,1), Hp+1, 1); repmat(umax, Hc, 1)]; % Limites superiores de g para X e u 
         
            
            %% 3.2. Definição dos Limites Operacionais
            La_max = 19;       % Limite máximo para o Nível do Anular (m).
            qm_max = 3.8;      % Limite máximo para a Vazão na coluna (m³/s).
            La_min = 5;        % Limite mínimo para o Nível do Anular (m).
            qm_min = 0.7;      % Limite mínimo para a Vazão na coluna (m³/s).
            H_max = 180;       % Limite superior de Head da bomba.
            H_min = 49;        % Limite inferior de Head da bomba.
            
            % Limites de Head baseados na vazão.
            [Hut,Hdt] = limites_DT_UP(obj.Envelope.limites,DadosProcesso(3)); % DadosProcesso(3) é qm
            ymax = [SP(1) * 1.05; Hut]; % Faixa superior para Nível (SP * 1.08) e Head (Hut).
            ymin = [SP(1) * 0.95; Hdt]; % Faixa inferior para Nível (SP * 0.92) e Head (Hdt)
           
            % Limites superiores e inferiores para as variáveis de decisão do solver (x0, du, ysp).
            % Estes são os limites 'lbx' e 'ubx' do solver.
            LimitesMax = [repmat([La_max; H_max; qm_max], Hp+1, 1); repmat(dumax, Hc, 1);  ymax(1); ymax(2)];
            LimitesMin = [repmat([La_min; H_min; qm_min], Hp+1, 1); repmat(-dumax, Hc, 1); ymin(1); ymin(2)];        
            
            %% 3.3. Preparação de Dados para o Solver
            % Busca a predição do estado no horizonte mais distante (Hp) para cálculo de erro.
            PredicaoHorizonteHp = obj.x0(nx*Hp + 1 : nx*(Hp+1));
            % Calcula o Head predito com base nas variáveis preditas em Hp.
            PredicaoHorizonteHead = PredicaoHorizonteHp(2);
            % Cálculo do erro (La, Head) entre o processo e a predição no horizonte Hp.
            erro_feedback = ([DadosProcesso(1);DadosProcesso(2)] - [PredicaoHorizonteHp(1);PredicaoHorizonteHead]);
            % Valor padrão para indicar que o solver ainda não foi executado 
            Feasible = 0.5; 
            % Vetor de parâmetros para o solver CasADi (input 'p').
            ParSolver = [DadosProcesso(1:3);obj.uk;erro_feedback;DadosProcesso(4:5);uRTO;obj.ModeloPreditor.data.a0];
                       
            %% 3.4. Lógica de Controle: Operação Manual vs. NMPC
            if t < obj.tempo_manual
                % Fase de operação manual: O processo se estabiliza.
                disp('Aguardando estabilizar o processo (operação manual)...');
                obj.uk = [obj.f_manual(t+1); obj.Zc_manual(t+1)]; % Usa valores do perfil manual.
            else
                % Fase NMPC: Executa o controlador preditivo.
                obj.contador = obj.contador + 1; % Incrementa o contador de passos do MPC.
                if obj.contador == PassoMPC 
                    tStar = tic; % Inicia a contagem do tempo de computação do solver.
                    sol_NMPC = obj.casadi_solver('x0', obj.x0, ...              % Palpite inicial
                                                'lbx', LimitesMin, ...          % Limites inferiores das variáveis de otimização
                                                'ubx', LimitesMax, ...          % Limites superiores das variáveis de otimização
                                                'lbg', ManipuladasLowLimit, ... % Limites inferiores das restrições 'g'
                                                'ubg', ManipuladasHighLimit, ...% Limites superiores das restrições 'g'
                                                'p', ParSolver);                % Parâmetros para o problema
                    Feasible = obj.casadi_solver.stats.success;
                    obj.x0 = full(sol_NMPC.x); % Atualiza o palpite inicial com a solução ótima encontrada
                    if Feasible<0.5
                       duk = obj.duk_anterior;        
                    else
                       duk = obj.x0(nx*(Hp+1)+1:nx*(Hp+1)+Hc*nu); % Incrementos de controle ótimos (Deltau_k).
                    end
                    obj.duk_anterior = duk(1:nu);     % Guardar a solução anterior
                    obj.uk      = obj.uk + duk(1:nu); % Nova ação de controle (ação anterior + primeiro incremento).
                    obj.fob     = full(sol_NMPC.f);   % Valor da Função Objetivo.
                    obj.ysp     = obj.x0(nx*(Hp+1)+Hc*nu+1:end); % Novos setpoints (se otimizados).
                    obj.Tsolver = toc(tStar);         % Tempo de computação do solver.
                    obj.contador = 0;                 % Reseta o contador para o próximo passo do MPC.
                end
            end
            
            %% 3.5. Atualização do estado da ESN após a decisão de controle
            ukk = feature_scaling(vertcat(obj.uk,DadosProcesso(4:5),DadosProcesso(1:3)),u_max(1:7),u_min(1:7));
            obj.ModeloPreditor.data.a0 = esquenta_ESNx(obj.ModeloPreditor.data,ukk,1);
            
            %% 3.6. Saída do Bloco NMPC
            % Organiza todas as informações relevantes em um único vetor de saída 'u'.
            u = [obj.uk;                                     % Ação de controle aplicada [nu x 1]
                        obj.ysp;                             % Setpoints atuais/ótimos [ny x 1]
                        ymax;                                % Limites superiores das saídas [ny x 1]
                        ymin;                                % Limites inferiores das saídas [ny x 1]
                        [La_max; H_max; qm_max; H_max];      % Limites máximos dos estados (La,H,qm) e H_max novamente? [nx+1 x 1] (se H_max já está em RestricoesMax_estados, redundante)
                        [La_min; H_min; qm_min; H_min];      % Limites mínimos dos estados (La,H,qm) e H_min novamente? [nx+1 x 1]
                        ManipuladasHighLimit(end-1:end);     % Limites superiores das variáveis manipuladas umax_op [nu x 1]
                        ManipuladasLowLimit(end-1:end);      % Limites inferiores das variáveis manipuladas umin_op [nu x 1]
                        obj.fob;                             % Valor da função objetivo [1x1]
                        obj.Tsolver;                         % Tempo de cálculo do solver [1x1]
                        Feasible;                            % Status de viabilidade [1x1]
                        PredicaoHorizonteHp;                 % Predição dos estados no final do horizonte anterior [nx x 1]
                        PredicaoHorizonteHead];              % Predição do Head no final do horizonte anterior [1x1]
            
        end
    end
end