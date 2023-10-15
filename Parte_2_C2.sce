// Sistema de EDOs Não-Linear

function dxdt = h(t,x)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    dxdt = [x3;
    x4;
    (33/9*sin(x1)+sin(x1-2*x2)+sin(2*x1-2*x2)*x3^2+8/9*sin(x1-x2)*x4^2)/(-3+cos(2*x1-2*x2));
    (-5.25*sin(2*x1-x2)+3.75*sin(x2)-9*sin(x1-x2)*x3^2-sin(2*x1-2*x2)*x4^2)/(-3+cos(2*x1-2*x2))];
endfunction

// Condições Iniciais (C2)

x10 = %pi/2;
x20 = %pi/10;
x30 = 0;
x40 = 0;

x0 = [x10;x20;x30;x40];

// Tempo de Simulação (C2)

t0 = 0; 
tf = 48;
res = 3000;

tempo = linspace(t0,tf,res);

// Solução do Sistema Não-Linear (M1)

SN_x = ode("rk",x0,t0,tempo,h);

// Gráfico da Energia Mecânica (M1) "Pensar em Variação!"

show_window(0)

g = 9.81;

theta1N = SN_x(1,:);
theta2N = SN_x(2,:);
theta1Np = SN_x(3,:);
theta2Np = SN_x(4,:);

i=1;
EM=[];
EM0=g^3*(2/9*theta1Np(i)*theta2Np(i)*cos(theta1N(i)-theta2N(i))+0.5*theta1Np(i)^2-21/18*cos(theta1N(i))+4/81*theta2Np(i)^2-2/9*cos(theta2N(i)));

while i<=res,
    EM($+1) = g^3*(2/9*theta1Np(i)*theta2Np(i)*cos(theta1N(i)-theta2N(i))+0.5*theta1Np(i)^2-21/18*cos(theta1N(i))+4/81*theta2Np(i)^2-2/9*cos(theta2N(i)))-EM0;
    i=i+1; end

plot2d(tempo, EM, 2)
title ("Energia Mecânica em Função do Tempo - M1")
xlabel("t (s)")
ylabel("EM (J.m/kg)")

show_window(1)

plot2d(tempo, EM, 2)
xlabel("t (s)")
ylabel("EM (J.m/kg)")
legends(leg="M1", 2,opt=2)

// Solução do Sistema Não-Linear (M2)

SN_x = ode("adams",x0,t0,tempo,h);

// Gráfico da Energia Mecânica (M2)

theta1N = SN_x(1,:);
theta2N = SN_x(2,:);
theta1Np = SN_x(3,:);
theta2Np = SN_x(4,:);

i=1;
EM=[];
EM0=g^3*(2/9*theta1Np(i)*theta2Np(i)*cos(theta1N(i)-theta2N(i))+0.5*theta1Np(i)^2-21/18*cos(theta1N(i))+4/81*theta2Np(i)^2-2/9*cos(theta2N(i)));

while i<=res,
    EM($+1) = g^3*(2/9*theta1Np(i)*theta2Np(i)*cos(theta1N(i)-theta2N(i))+0.5*theta1Np(i)^2-21/18*cos(theta1N(i))+4/81*theta2Np(i)^2-2/9*cos(theta2N(i)))-EM0;
    i=i+1; end

plot2d(tempo, EM, 5)
title ("Energia Mecânica em Função do Tempo - M1 & M2")
xlabel("t (s)")
ylabel("EM (J.m/kg)")
legends(leg="M2", 5,opt=3)

show_window(2)

plot2d(tempo, EM, 5)
title ("Energia Mecânica em Função do Tempo - M2")
xlabel("t (s)")
ylabel("EM (J.m/kg)")


