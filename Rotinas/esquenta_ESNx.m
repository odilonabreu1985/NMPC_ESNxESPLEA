function proximo_estado = esquenta_ESNx(modelo_esn, entradas_normalizadas,QtePassosESN)
leakrate = modelo_esn.gama;
estado_reservatorio_esn = modelo_esn.a0;
for i=1:QtePassosESN
    % calcula novos estados da ESN
    xk_1 = modelo_esn.Wrr*estado_reservatorio_esn +  modelo_esn.Wir*entradas_normalizadas + modelo_esn.Wbr; 
    estado_reservatorio_esn = (1-leakrate)*estado_reservatorio_esn + leakrate*tanh(xk_1);
end
proximo_estado = estado_reservatorio_esn; 
