function [sys,x0]=envelope_LEA(t,x,u,flag,Par,qm_x0,head_x0);
% Animação para atualização do ponto de operação em tempo real
global PontoOp
global AlvoMapa
global TituloVazao
global VazaoAlvo
global HeadAlvo
global PontoOp1

%==================================================================================
% Inicialização
if flag==0
    % Inicializa para S-Function
    ninput=3;
    nout=0;
    x0=[];
    sys= [0;size(x0,1); nout; ninput;0;0];
    monta_mapa_LEA(Par);
    VazaoAlvo = 3.56;
    HeadAlvo = 117;
    PontoOp=plot(qm_x0*3600,head_x0,'ko','MarkerSize',8);
    scatter(VazaoAlvo,HeadAlvo,'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 1 0],...
        'LineWidth',3)
    %AlvoMapa=plot(VazaoAlvo,HeadAlvo,'x','MarkerSize',10,'Color',[0 1 0]);    % Plota alvo incial no mapa
    %TituloVazao=title(strcat("VazaoOleoEstimada=",num2str(qm_x0,'%.1f')," (*) SetPoint=",num2str(head_x0,'%.1f')," Erro = ",num2str(qm_x0-VazaoAlvo,'%.3f')," [m³/h]"));
    TituloVazao=title(strcat("H_{sp}*=",num2str(head_x0,'%.1f'),"; qm=",num2str(qm_x0,'%.3f'),"; ErroTarget = ",num2str(qm_x0-VazaoAlvo,'%.3f')));
    PontoOp1=plot(qm_x0*3600,head_x0,'k*','MarkerSize',8);
    
    %==================================================================================
elseif flag==2
    qm1=  PontoOp.XData;   % Ponto da vazao, antes de atualizar
    Head1=PontoOp.YData;   % Ponto de operacao Head, antes de atualizar
    qm=u(1);               % Vazao na coluna de producao estimada (m³/h)
    Head=u(2);             % Head estimado (m)
    HeadOtimo = u(3);
    plot([qm1  qm],[Head1  Head],'k','MarkerSize',12)
    set(PontoOp, 'xdata', qm, 'ydata',Head);
    set(PontoOp1,'xdata', qm, 'ydata',HeadOtimo);
    TituloVazao.String=(strcat("H_{sp}*=",num2str(HeadOtimo,'%.1f'),"; qm=",num2str(qm,'%.3f'),"; ErroTarget = ",num2str(qm-VazaoAlvo,'%.3f')));
    drawnow;
    sys=[];
end
%==================================================================================
%==================================================================================
%==================================================================================
function monta_mapa_LEA(Par)
close all
% =========================================================================
% Correções de viscosidade Turzo
% =========================================================================
H_aguabep = 4.330800000000000e+02; %ft
Q_aguabep = 4.401900000000000e+02; %bpd
y = -112.1374+6.6504*log(H_aguabep)+12.8429*log(Q_aguabep);
Q = (exp((39.5276+26.5605*log(Par.visc*1000)-y)/51.6565)); %Pa.s to Cp;
Cq = (1.0-4.0327e-03*Q-1.7240e-04*Q^2);
CH = (1.0-4.4723e-03*Q) -(4.1800e-05*Q^2); % 80%
Qdown0 =  Par.q0_min;
Qup0 =    Par.q0_max;
H_down_0 = Par.Head(1)*Qdown0^4 +Par.Head(2)*Qdown0^3 + Par.Head(3)*Qdown0^2 + Par.Head(4)*Qdown0 + Par.Head(5);
H_up_0 =   Par.Head(1)*Qup0^4 +  Par.Head(2)*Qup0^3 + Par.Head(3)*Qup0^2 + Par.Head(4)*Qup0 + Par.Head(5);
Frequencias = 25:(Par.fMax-25)/2000:Par.fMax;

H_down = CH.*H_down_0.*(Frequencias./Par.f0).^2;
H_up   = CH.*H_up_0.*(Frequencias./Par.f0).^2;
Qdown30 = Qdown0./Cq.*Frequencias./Par.f0;
Qup30 =   Qup0./Cq.*Frequencias./Par.f0;

Q_max_graParh = 5/3600;
vazoes = (0:Q_max_graParh/200:Q_max_graParh)./Cq;
Head_60Hz = CH.*(Par.Head(1).*vazoes.^4 +  Par.Head(2).*vazoes.^3 + Par.Head(3).*vazoes.^2 + Par.Head(4).*vazoes + Par.Head(5));
Head_30Hz = Head_60Hz.*(30/Par.f0).^2;
Q_30Hz = vazoes.*(30/60);
Head_35Hz = Head_60Hz.*(35/Par.f0).^2;
Q_35Hz = vazoes.*(35/60);
Head_40Hz = Head_60Hz.*(40/Par.f0).^2;
Q_40Hz = vazoes.*(40/60);
Head_45Hz = Head_60Hz.*(45/Par.f0).^2;
Q_45Hz = vazoes.*(45/60);
Head_50Hz = Head_60Hz.*(50/Par.f0).^2;
Q_50Hz = vazoes.*(50/60);
Head_55Hz = Head_60Hz.*(55/Par.f0).^2;
Q_55Hz = vazoes.*(55/60);

% Definição da região de incerteza
Qdown0_inc = Par.q0_min.*1.10; % Vazão mínima
Qup0_inc =   Par.q0_max.*0.90; % Vazão máxima
H_down_0_inc = Par.Head(1).*Qdown0_inc.^4 +  Par.Head(2).*Qdown0_inc.^3 + Par.Head(3).*Qdown0_inc.^2 + Par.Head(4).*Qdown0_inc + Par.Head(5);
H_up_0_inc   = Par.Head(1).*Qup0_inc.^4   +  Par.Head(2).*Qup0_inc.^3 + Par.Head(3).*Qup0_inc.^2 + Par.Head(4).*Qup0_inc + Par.Head(5);
% Corrige o Head com a lei de afinidades e a frequências
H_down_inc = CH.*H_down_0_inc.*(Frequencias./Par.f0).^2;
H_up_inc   = CH.*H_up_0_inc.*(Frequencias./Par.f0).^2;
% Corrige as vazões com a lei de afinidades e a frequência
Qdown30_inc = Qdown0_inc./Cq.*Frequencias./Par.f0;
Qup30_inc =   Qup0_inc./Cq.*Frequencias./Par.f0;


% =========================================================================
% Definição dos Polígonos de UPthrust e Downthrust
% =========================================================================
[~,~,Aii] = polyxpoly(Q_30Hz,Head_30Hz,Qdown30,H_down);
[~,~,Bii] = polyxpoly(Q_30Hz,Head_30Hz,Qup30,H_up);
[~,~,Cii] = polyxpoly(vazoes,Head_60Hz,Qup30,H_up);
[Dx,~,Dii] = polyxpoly(vazoes,Head_60Hz,Qdown30,H_down);
% Definição dos polígonos de incerteza
[~,~,Eii] = polyxpoly(Q_30Hz,Head_30Hz,Qdown30_inc,H_down_inc);
[~,~,Fii] = polyxpoly(Q_30Hz,Head_30Hz,Qup30_inc,H_up_inc);
[~,~,Gii] = polyxpoly(vazoes,Head_60Hz,Qup30_inc,H_up_inc);
[~,~,Hii] = polyxpoly(vazoes,Head_60Hz,Qdown30_inc,H_down_inc);

Vazao_down_60Hz = (0:Dx/99:Dx);
Head_60Hz_down = CH.*(Par.Head(1).*(Vazao_down_60Hz).^4 + Par.Head(2).*(Vazao_down_60Hz).^3 + Par.Head(3).*(Vazao_down_60Hz).^2 + Par.Head(4).*(Vazao_down_60Hz) + Par.Head(5))*(60/Par.f0).^2;
% =========================================================================
% Definição Upthrust e Downthrust
% =========================================================================
% Polígono de Downtrust
PP1 =     [[Q_30Hz(1:Aii)*3600; Head_30Hz(1:Aii)]'; [(Qdown30(Aii(2):Dii(2)-1)*3600)'     H_down(Aii(2):Dii(2)-1)'];    [Vazao_down_60Hz(end:-1:2)*3600 ; Head_60Hz_down(end:-1:2)]'];
PP1_inc = [[Q_30Hz(1:Eii)*3600; Head_30Hz(1:Eii)]'; [(Qdown30_inc(Eii(2):Hii(2)-1)*3600)' H_down_inc(Eii(2):Hii(2)-1)'];[Vazao_down_60Hz(end:-1:2)*3600 ; Head_60Hz_down(end:-1:2)]'];

% Polígono de UPthrust
PP2 = [[Q_30Hz(end)*3600; Head_30Hz(end)]';[Q_35Hz(end)*3600 ; Head_35Hz(end)]' ; ...
    [Q_40Hz(end)*3600 ; Head_40Hz(end)]';[Q_45Hz(end)*3600 ; Head_45Hz(end)]'; ...
    [Q_50Hz(end)*3600 ; Head_50Hz(end)]';[Q_55Hz(end)*3600 ; Head_55Hz(end)]';...
    [vazoes(end:-1:Cii(1))*3600 ; Head_60Hz(end:-1:Cii(1))]';...
    [(Qup30(Cii(2):-1:Bii(2))*3600)' H_up(Cii(2):-1:Bii(2))'];...
    [(Q_30Hz(Bii:end)*3600)' (Head_30Hz(Bii:end))']];

PP2_inc = [[Q_30Hz(end)*3600; Head_30Hz(end)]';[Q_35Hz(end)*3600 ; Head_35Hz(end)]' ; ...
    [Q_40Hz(end)*3600 ; Head_40Hz(end)]';[Q_45Hz(end)*3600 ; Head_45Hz(end)]'; ...
    [Q_50Hz(end)*3600 ; Head_50Hz(end)]';[Q_55Hz(end)*3600 ; Head_55Hz(end)]'; ...
    [vazoes(end:-1:Gii(1))*3600 ; Head_60Hz(end:-1:Gii(1))]'; ...
    [(Qup30_inc(Gii(2):-1:Fii(2))*3600)' H_up_inc(Gii(2):-1:Fii(2))'];[(Q_30Hz(Fii:end)*3600)' (Head_30Hz(Fii:end))']];

pgon1 = polyshape(PP1,'Simplify',false);
pgon1_inc = polyshape(PP1_inc,'Simplify',false);
pgon2 = polyshape(PP2,'Simplify',false);
pgon2_inc = polyshape(PP2_inc,'Simplify',false);

figure1 = figure('position',[1    300   689   451],'color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
p1 = plot(vazoes*3600,Head_60Hz,'k-','LineWidth',2);
hold on
p2 = plot(Q_30Hz*3600,Head_30Hz,'k-','LineWidth',2);
p2 = plot(Q_35Hz*3600,Head_35Hz,'--','LineWidth',2,'color',0.8*[1 1 1]);
p2 = plot(Q_40Hz*3600,Head_40Hz,'--','LineWidth',2,'color',0.65*[1 1 1]);
p2 = plot(Q_45Hz*3600,Head_45Hz,'--','LineWidth',2,'color',0.8*[1 1 1]);
p2 = plot(Q_50Hz*3600,Head_50Hz,'--','LineWidth',2,'color',0.65*[1 1 1]);
p2 = plot(Q_55Hz*3600,Head_55Hz,'--','LineWidth',2,'color',0.8*[1 1 1]);
p3 = plot(pgon1,'FaceColor','black','FaceAlpha',0.1,'LineWidth',2);
p3 = plot(pgon1_inc,'FaceColor','black','FaceAlpha',0.1,'LineStyle','none'); % Polígono de Incertaza
p4 = plot(pgon2,'FaceColor','black','FaceAlpha',0.1,'LineWidth',2);
p4 = plot(pgon2_inc,'FaceColor','black','FaceAlpha',0.1,'LineStyle','none'); % Polígono de Incertaza
text(vazoes(1)*3600-0.3,Head_60Hz(1),'60Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_55Hz(1),'55Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_50Hz(1),'50Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_45Hz(1),'45Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_40Hz(1),'40Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_35Hz(1),'35Hz','color','black','FontSize',14)
text(vazoes(1)*3600-0.3,Head_30Hz(1),'30Hz','color','black','FontSize',14)
txt = {'Downthrust'};
ttt = text(0.151132638810807,80.33563466291244,0,txt,'FontSize',20,'Color',[0 0 0]);
set(ttt,'Rotation',85);
txt = {'Upthrust'};
ttt = text(3.099869957565817,40.531533665598232,0,txt,'FontSize',20,'Color',[0 0 0]);
set(ttt,'Rotation',40);
ylim([0 200])
xlim([-0.9 5.5])
xlabel({'q_{m}/(m³h^{-1})'},'FontWeight','bold','FontAngle','italic','FontSize',14)
ylabel({'H/(m)'},'FontWeight','bold','FontAngle','italic','FontSize',14);
box('on')
set(axes1,'FontSize',16,'FontWeight','bold');
set(figure1,'Units','Inches');







  

