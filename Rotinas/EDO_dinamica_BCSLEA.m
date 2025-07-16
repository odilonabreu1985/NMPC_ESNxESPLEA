function rhs=EDO_dinamica_BCSLEA(xmk,ukk,Par)
% =========================================================================
% Define as entradas 
% =========================================================================
    freq =   ukk(1);
    valv =   ukk(2);
    pman =   ukk(3);
    preser = ukk(4);
% =========================================================================
% Define os estados 
% =========================================================================
    La =     xmk(1);
    pchoke = xmk(2);
    qm =     xmk(3);
    visc =   xmk(4);
    IP =     xmk(5);
% =========================================================================
% Variáveis algébricas  
% =========================================================================
    pfundo = Par.Vol_1.rho_1*Par.g*(Par.Vol_1.h_r-Par.Vol_1.h_p)+Par.Vol_1.rho_1*Par.g*La;
    qc_choke = ((valv).*((pchoke-pman).^0.5))*Par.k_choke;
    qreser = (IP*(preser - pfundo));
% =========================================================================
% Define a dinâmica da válvula choke, Fricção e Deltap da bomba, Pressao de
% fundo e Vazão de reservátorio
% =========================================================================
    Re_1 = (4*Par.Vol_1.rho_1*qm/1000)/(0.219*pi*visc);
    Re_2 = (4*Par.Vol_2.rho_2*qm/1000)/(0.219*pi*visc);

    F_1 = ((64/Re_1)*qm^2/1000*Par.Vol_1.rho_1)/(2*pi*Par.Vol_1.r1^3);
    F_2 = ((64/Re_2)*qm^2/1000*Par.Vol_2.rho_2)/(2*pi*Par.Vol_2.r2^3);
    
    H_aguabep = 4.330800000000000e+02; %ft
    Q_aguabep = 4.401900000000000e+02; %bpd
    y = -112.1374+6.6504*log(H_aguabep)+12.8429*log(Q_aguabep);
    Q = (exp((39.5276+26.5605*log(visc*1000)-y)/51.6565)); %Pa.s to Cp;
    Cq = (1.0-4.0327e-03*Q-1.7240e-04*Q^2);
    CH = (1.0-4.4723e-03*Q) -(4.1800e-05*Q^2); % 80%
    q_0a = qm/Cq*(Par.f0/freq);
    H_0 =  Par.Head(1)*q_0a^4 +  Par.Head(2)*q_0a^3 + Par.Head(3)*q_0a^2 + Par.Head(4)*q_0a + Par.Head(5);
    H = (CH*H_0*(freq/Par.f0)^2);
    dp = H*Par.Vol_1.rho_1*Par.g; 
% =========================================================================
% Resolve a E.D.O do sistema BCS 
% =========================================================================
    f1 = (qreser - qm)/Par.Anular;
    f2 = (Par.Vol_2.B2/Par.Vol_2.V2)*(qm-qc_choke);
    f3 = (Par.A_bar/(Par.rho_bar*Par.l_bar))*(pfundo - pchoke - F_1 - F_2 - ...
         (Par.Vol_1.rho_1*Par.g*(Par.Vol_1.h_r-Par.Vol_1.h_p))-(Par.Vol_2.rho_2*Par.g*(Par.Vol_1.h_p-Par.Vol_2.h_c))+ (dp));
    f4 = 0;
    f5 = 0;
    rhs = [f1;f2;f3;f4;f5];
   
        
