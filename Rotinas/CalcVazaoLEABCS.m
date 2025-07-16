function  [sys, x0]  = CalcVazaoLEABCS(t,x,u,flag,DeltaT,QIni)
% Função para calcular a produção de óleo acumulada ao longo da simulação
%---------------------------------------------------------------------------------------------
% flag = 0 --> Define condições iniciais
if flag == 0
   ninput=1;			  % Num de entradas = Freq e PChegada
   nout=1;				   % Num. de saida = Vazão, Produção de Óleo e Produção de Óleo Acumulada
   x0 = QIni;		% Inicializa vetor de estados (coluna), no caso, Vazão Inicial, Produção=0 e Produção Acumulada=0
   sys= [0;size(x0,1); nout; ninput;0;0];  % Depois que defini, nunca mudei desde 199x! :-)
%--------------------------------------------------------------------------------------------
% flag = 2 --> Retorna estados do sistema
elseif flag == 2
    Vazao=u(1)/3600;                                % Resgata valor da Frequência na entrada do bloco
    ProducaoAcumulada=x(1);                         % Resgata valor da produção acumulada
    Producao=Vazao;                                 % Vazão Instantânea em m3/s
    Producao=Producao*DeltaT;                       % Produção (m3) na janela de tempo (Ts padrão = 10s)
    ProducaoAcumulada=ProducaoAcumulada+Producao;   % Incrementa produção acumulada com a produção atual
    sys=ProducaoAcumulada;                  	    % Retorna o valor total da produção e da produção acumulada
    %--------------------------------------------------------------------------------------------
% flag = 3 --> Retorna vetor de saida
elseif flag == 3
	sys= x;				% A saída é o próprio vetor de estados com a Vazão, Produção e a Produção Acumulada
end
%--------------------------------------------------------------------------------------------


