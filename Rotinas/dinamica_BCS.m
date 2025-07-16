    import casadi.*
% =========================================================================
% Define as variáveis simbólicas nos estados e nos parâmetros   
% =========================================================================
    x = MX.sym('x',3);   
    u = MX.sym('u',5);
% =========================================================================
% Define as entradas 
% =========================================================================
    freq = u(1);
    valv = u(2);
    pman = u(3);
    preser = u(4);
    Temp_fluido = u(5);

% =========================================================================
% Define os estados 
% =========================================================================
    La =      x(1);
    pchoke =  x(2);
    qm =      x(3);

% =========================================================================
% Variaveis calculadas 
% =========================================================================
    visc = (30.242*exp(-0.024*Temp_fluido))/1000;
    pressao_fundo = Par.Vol_1.rho_1*Par.g*(Par.Vol_1.h_r-Par.Vol_1.h_p)+Par.Vol_1.rho_1*Par.g*La;
    qc_choke = ((valv).*((pchoke-pman).^0.5))*Par.k_choke;
    qc_reser = (Par.IP*(preser - pressao_fundo));
    pin =      La*Par.Vol_1.rho_1*Par.g;                                                       

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
   %CH = (1.0-3.6800e-03*Q) -(4.3600e-05*Q^2); % 60%
    CH = (1.0-4.4723e-03*Q) -(4.1800e-05*Q^2); % 80%
   %CH = (1.0-7.007e-03*Q) -(1.4100e-05*Q^2); % 100%
   %CH = (1.0-9.0100e-03*Q) +(1.3100e-05*Q^2); % 120%
    q_0a = qm/Cq*(Par.f0/freq);
   %H_0 =  Par.Head(1)*q_0a^2 +  Par.Head(2)*q_0a +  Par.Head(3);
    H_0 =  Par.Head(1)*q_0a^4 +  Par.Head(2)*q_0a^3 + Par.Head(3)*q_0a^2 + Par.Head(4)*q_0a + Par.Head(5);
    H = (CH*H_0*(freq/Par.f0)^2);
    dp = H*Par.Vol_1.rho_1*Par.g; 
    P_0_aux = Par.Pot(1)*q_0a^3 + Par.Pot(2)*q_0a^2 + Par.Pot(3)*q_0a + Par.Pot(4);
    BHP = P_0_aux*(freq/Par.f0)^3;               % Watts
    P_hp = BHP*0.001341022;                      % HP
    Pc = (Par.Vol_1.rho_1*H*q_0a);               % Eixo da bomba (potencia hidráulica w)
    rend = Par.Eff(1)*qc_choke^4+ Par.Eff(2)*qc_choke^3+ Par.Eff(3)*qc_choke.^2 + Par.Eff(4)*qc_choke + Par.Eff(5);
    rend2 = Pc/BHP;
    I = (Par.I_np/Par.P_np)*BHP;                         % I_np [A] e P_np [W]
    WC = 0/100;          %BSW
    Co = 69.82*5.10;                                     %($/bpd) to (R$/bpd)  
    Ce = 60;                                             %(R$/kWh) preço da potência eletrica (???) 
   %L = -Co*((1-WC)*(qc_choke*543439.7)) + Ce*BHP*0.001; % maximizar producao e minimizar potência (Lucro = receita-custo) 
    L = -Co*((1-WC)*(qc_choke*543439.7));                % maximizar a produção
   %L =  Ce*BHP;                                         % minimizar a potência 
   %L =  -rend;                                          % maximizar a rendimento

% =========================================================================
% Define a EDO 
% =========================================================================
    f1 = (qc_reser - qm)/Par.Anular;
    f2 = (Par.Vol_2.B2/Par.Vol_2.V2)*(qm-qc_choke);
    f3 = (Par.A_bar/(Par.rho_bar*Par.l_bar))*(pressao_fundo - pchoke - F_1 - F_2 - ...
         (Par.Vol_1.rho_1*Par.g*(Par.Vol_1.h_r-Par.Vol_1.h_p))-(Par.Vol_2.rho_2*Par.g*(Par.Vol_1.h_p-Par.Vol_2.h_c))+ (dp));
    rhs = [f1;f2;f3];
% =========================================================================
% Define estados e saídas do sistema BCS 
% =========================================================================
    f_estados = Function('f_estados', {x,u},{rhs},{'x','u'},{'rhs'});
    y = [pin;H;P_hp;qc_choke;qc_reser;visc;I;L;rend;rend2]; 
    f_saida = Function('f_saida',{x,u},{y,pin,H,P_hp,qc_choke,qc_reser,visc,I,L,rend,rend2},{'x','u'},{'y','pin','H','P_hp','qc_choke','qc_reservatorio','visc','I','L','rend','rend2'}); % Sistema de Eq. Algebricas variaveis de saída
    BCS.estados = f_estados(x,u);
    BCS.saida = f_saida;
% =========================================================================
% Define estados e saídas do sistema BCS 
% =========================================================================
    f_estados = Function('f_estados', {x,u},{rhs},{'x','u'},{'rhs'});
    y = [pin;H;dp;P_hp;qc_choke;qc_reser;visc;Par.IP]; 
    f_saida = Function('f_saida',{x,u},{y,pin,H,dp,P_hp,qc_choke,qc_reser,visc,Par.IP},{'x','u'},{'y','pin','H','dp','P_hp','qc_choke','qc_reser','visc','Par.IP'}); % Sistema de Eq. Algebricas variaveis de saída
    BCS.estados = f_estados(x,u);
    BCS.saida = f_saida;
    ode = struct;     
    ode.x   = x;   
    ode.p   = u;   
    ode.ode = BCS.estados; %dinamica_BCS(states,inputs,Par);
    options = struct('tf',Par.Ts,'t0',0);
    intg = integrator('intg','cvodes',ode, options);
    res = intg('x0',x,'p',u);
    Eq_medicao = Function('Eq_medicao',{x,u},{y([1,2,3,4,5,6,7,8])},{'states','inputs'},{'y'});
    Sim_BCS_LEA = Function('Sim_BCS_LEA',{x,u},{res.xf,Eq_medicao(res.xf,u)},{'xk_1','uk_1'},{'xk','yk'});

    