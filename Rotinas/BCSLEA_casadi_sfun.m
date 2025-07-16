function [sys,x0] = BCSLEA_casadi_sfun(t,xend,u,flag,xmk,Sim_BCS_LEA)
% =========================================================================
% flag = 0 --> Define condições iniciais do processo   
% =========================================================================

if flag == 0
    nu = 5;             % número das entradas do processo                                    
    nx = size(xmk,1);   % número dos estados  do processo
	ninput=nu;			
	nout=nx;					
	x0 = xmk;	  % define as condições iniciais da ESN
    sys= [0;size(x0,1);nout;ninput;0;0];
% =========================================================================
% flag = 2 --> Retorna os estados do sistema 
% =========================================================================
elseif abs(flag) == 2
    [xk,ypk]=Sim_BCS_LEA(xend(1:3),u);
 	sys =[full(xk);full(ypk)]';
% =========================================================================
% flag = 3 --> Retorna vetor de saida
% =========================================================================

elseif flag == 3
	sys= xend;				  % Retorna apenas o novo y
end
  