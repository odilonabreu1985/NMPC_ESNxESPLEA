function  [sys, x0]  = CalcVazaoLEABCS(t,x,u,flag,DeltaT,QIni)
% Fun��o para calcular a produ��o de �leo acumulada ao longo da simula��o
%---------------------------------------------------------------------------------------------
% flag = 0 --> Define condi��es iniciais
if flag == 0
   ninput=1;			  % Num de entradas = Freq e PChegada
   nout=1;				   % Num. de saida = Vaz�o, Produ��o de �leo e Produ��o de �leo Acumulada
   x0 = QIni;		% Inicializa vetor de estados (coluna), no caso, Vaz�o Inicial, Produ��o=0 e Produ��o Acumulada=0
   sys= [0;size(x0,1); nout; ninput;0;0];  % Depois que defini, nunca mudei desde 199x! :-)
%--------------------------------------------------------------------------------------------
% flag = 2 --> Retorna estados do sistema
elseif flag == 2
    Vazao=u(1)/3600;                                % Resgata valor da Frequ�ncia na entrada do bloco
    ProducaoAcumulada=x(1);                         % Resgata valor da produ��o acumulada
    Producao=Vazao;                                 % Vaz�o Instant�nea em m3/s
    Producao=Producao*DeltaT;                       % Produ��o (m3) na janela de tempo (Ts padr�o = 10s)
    ProducaoAcumulada=ProducaoAcumulada+Producao;   % Incrementa produ��o acumulada com a produ��o atual
    sys=ProducaoAcumulada;                  	    % Retorna o valor total da produ��o e da produ��o acumulada
    %--------------------------------------------------------------------------------------------
% flag = 3 --> Retorna vetor de saida
elseif flag == 3
	sys= x;				% A sa�da � o pr�prio vetor de estados com a Vaz�o, Produ��o e a Produ��o Acumulada
end
%--------------------------------------------------------------------------------------------


