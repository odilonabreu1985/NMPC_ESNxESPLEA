function [Hut,Hdt] = limites_DT_UP(limites,Q)
% Identifica se a vaz�o est� dentro da faixas min e m�x do HEAD
% F_sat_dt=max(min(limites.Down(1,:)),min(max(limites.Down(1,:)),Q));
% F_sat_ut=max(min(limites.Up(1,:))  ,min(max(limites.Up(1,:)),Q));
F_sat_dt = min(max(limites.Down(1,:)), max(min(limites.Down(1,:)),Q));
F_sat_ut = min(max(limites.Up(1,:))  , max(min(limites.Up(1,:)),Q));

% C�lcula interpola��o linear entre os pontos existente para a vaz�o de entrada
Hut=interp1(limites.Up(1,:),limites.Up(2,:),F_sat_ut,'linear');
Hdt=interp1(limites.Down(1,:),limites.Down(2,:),F_sat_dt,'linear');


end

