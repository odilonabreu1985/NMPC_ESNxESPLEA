%% ===========================%  Par�metros do sistema BCS LEA ==================================================
function Par = dados_LEA_Par(Ts)
    Par.Ts = Ts;
%Volume de controle 1
    Par.Vol_1.B1 = 1.800e+09;                                  % [Pa]     Compressibilidade
    Par.Vol_1.rho_1 = 836.8898*1.02;                                % [kg/m^3] Densidade do fluido obtida pela m�dia dos experimentos no coriolis          
    Par.Vol_1.r1 = 0.11;                                       % [m]      Raio do tubo 
    Par.Vol_1.l1 =9.3;                                         % [m]      Comprimento do reservatorio at� a intake
    Par.Vol_1.h_r = 32;                                        % [m]      Profundidade do reservat�rio
    Par.Vol_1.h_p = 22.7;                                      % [m]      Profundidade do conjunto BCS
    Par.Vol_1.V1 = Par.Vol_1.l1*pi*Par.Vol_1.r1^2;             % [m^3]    C�lculo do volume 
    
%Volume de controle 2
    Par.Vol_2.B2 = 1.800e+09;                                 % [Pa]     Compressibilidade
    Par.Vol_2.rho_2 = Par.Vol_1.rho_1;                        % [kg/m^3] Densidade do fluido obtida pela m�dia dos experimentos no coriolis          
    Par.Vol_2.r2 = 0.0375;                                     % [m]      Raio do tubo 
    Par.Vol_2.l2 = 22.7;                                       % [m]      Comprimento da intake at� a choke
    Par.Vol_2.h_c = 0;                                         % [m]      Profundidade do choke
    Par.Vol_2.V2 = Par.Vol_2.l2*pi*Par.Vol_2.r2^2;             % [m^3]    C�lculo do volume 
    
    % Outras
    Par.Pman =    1.7;    % [bar] Press�o de Manifold 1 atm
    Par.Resv =    2.388125181198120*0.57;%2.05;                                      % [Pa] M�ximo do valor experimental
    Par.k_choke = 2.92e-08; 
    Par.visc =    0.012;                                       % [Pa.s] Viscosidade din�mica
    Par.g =       9.81;                                        % [m/s^2] Gravidade
    Par.IP =      3.3e-08;                                    % [m^3/(s*Pa)]Indice de Produtividade

    Par.fMax = 60;                                             % [Hz]   Frequ�ncia m�xima admiss�vel
    Par.fMin = 30;                                             % [Hz]   Frequ�ncia m�nima admiss�vel
    Par.f0 =   60;                                             % Frequecia inicial
    Par.Anular =  0.033595;                                    % Area do anular LEA
    Par.I_np = 25;                                             % [A] corrente nominal  
    Par.P_np = 13422.5982;                                     % [W] ou 18 [hp] obtida do Q_bep
    Par.q0_min = 0.5/3600; % Vaz�o m�nima de refer�ncia
    Par.q0_max = 4/3600;   % Vaz�o m�xima de refer�ncia
    Par.Tintake = 42;
    % Valores m�dios
    Par.l_bar = (Par.Vol_1.l1 + Par.Vol_2.l2)/2;
    Par.r_bar = ((Par.Vol_1.r1*Par.Vol_1.l1)+(Par.Vol_2.r2*Par.Vol_2.l2))/(Par.Vol_1.l1+Par.Vol_2.l2);
    Par.A_bar =  pi*Par.r_bar^2;
    Par.rho_bar = ((Par.Vol_1.rho_1*Par.Vol_1.V1)+(Par.Vol_2.rho_2*Par.Vol_2.V2))/(Par.Vol_1.V1+Par.Vol_2.V2);

    % Valores admensionais
    Par.L_poco = Par.Vol_1.h_r;
    Par.DeltaP_max =15.534099888869710;
    Par.Delta_vazao = 4*1000/3600;
    
    %Curva da bomba
    Par.Pot  = [3.14;3.41;3.57;3.60]*745.7; % Hp to Watts
    Par.Eff =  [-38000037638726.0,76885738063.1052,-104037865.052405,92277.4570774092,-0.000102311045885796];
    Par.Head = [77666347302028.2;-234691030380.973;156228165.550398;-63532.2711035997;189.799864460968]; %corrigido fminsearch (22_09_2021)

    %C�lculos auxiliares head
    %head = (X_BCS(2:2,end)-(X_BCS(1:1,end)*Par.Vol_1.rho_1*Par.g)+(Par.Vol_2.rho_2*Par.g*(Par.Vol_1.h_p-Par.Vol_2.h_c))+(((64/(4*Par.Vol_2.rho_2*X_BCS(3:3,end)/1000)/(0.219*pi*Par.visc))*X_BCS(3:3,end)^2/1000*Par.Vol_2.rho_2)/(2*pi*Par.Vol_2.r2^3)))/(Par.Vol_2.rho_2*Par.g);


   