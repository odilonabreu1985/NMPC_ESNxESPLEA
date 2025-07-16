function [u_max,u_min,y_max,y_min] = normalizar_BCS_ESN()
    % --- Limites para Variáveis de Controle (Entradas u) ---
    lim_freq_min = 40;    % Hz
    lim_freq_max = 60;    % Hz
    lim_choke_min = 0;    % %
    lim_choke_max = 100;  % %
    
    % --- Limites para Parâmetros do Processo (podem ser distúrbios ou variáveis exógenas) ---
    lim_pman_min = 0.5;   % bar (ou unidade consistente) - Pressão do manifold
    lim_pman_max = 5;     % bar
    lim_pres_min = 1;     % bar - Pressão do reservatório
    lim_pres_max = 2.3;   % bar
    
    % --- Limites para Saídas do Processo (Variáveis y) / Estados Principais ---
    lim_La_min = 0;   % m (ou unidade consistente) - Nível no anular
    lim_La_max = 20;  % m
    lim_H_min = 28;   % m (ou unidade consistente) - Head da bomba
    lim_H_max = 200;  % m
    lim_qm_min = 0;   % m³/h (ou unidade consistente) - Vazão de produção
    lim_qm_max = 5;   % m³/h
    
    u_max = [lim_freq_max;lim_choke_max;lim_pman_max;lim_pres_max;lim_La_max; lim_H_max; lim_qm_max];
    u_min = [lim_freq_min;lim_choke_min;lim_pman_min;lim_pres_min;lim_La_min; lim_H_min; lim_qm_min];
    
    y_max = [lim_La_max; lim_H_max; lim_qm_max];
    y_min = [lim_La_min; lim_H_min; lim_qm_min];
end

