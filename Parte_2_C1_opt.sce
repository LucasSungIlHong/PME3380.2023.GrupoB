// Sistema de EDOs Linear

function dxdt = f(t,x)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    dxdt = [x3; x4; (-7/3)*x1+x2; (21/4)*x1-(9/2)*x2];
endfunction

// Condições Iniciais (C1)

x10 = 0.21;
x20 = 0.21;
x30 = 0;
x40 = 0;

x0 = [x10;x20;x30;x40];

// Tempo de Simulação (C1)

t0 = 0; 
tf = 18;
res = 1500;

tempo = linspace(t0,tf,res);

// Solução do Sistema Linear (M1)

S_x = ode("rk",x0,t0,tempo,f);

// Gráfico da Energia Mecânica (M1) "Pensar em Variação!"

show_window(0)

g = 9.81;

theta1 = S_x(1,:);
theta2 = S_x(2,:);
theta1p = S_x(3,:);
theta2p = S_x(4,:);

i=1;
EM=[];
EM0=g^3*(2/9*theta1p(i)*theta2p(i)*cos(theta1(i)-theta2(i))+0.5*theta1p(i)^2-21/18*cos(theta1(i))+4/81*theta2p(i)^2-2/9*cos(theta2(i)));

while i<=res,
    EM($+1) = g^3*(2/9*theta1p(i)*theta2p(i)*cos(theta1(i)-theta2(i))+0.5*theta1p(i)^2-21/18*cos(theta1(i))+4/81*theta2p(i)^2-2/9*cos(theta2(i)))-EM0;
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

S_x = ode("adams",x0,t0,tempo,f);

// Gráfico da Energia Mecânica (M2)

theta1 = S_x(1,:);
theta2 = S_x(2,:);
theta1p = S_x(3,:);
theta2p = S_x(4,:);

i=1;
EM_=[];
EM_0=g^3*(2/9*theta1p(i)*theta2p(i)*cos(theta1(i)-theta2(i))+0.5*theta1p(i)^2-21/18*cos(theta1(i))+4/81*theta2p(i)^2-2/9*cos(theta2(i)));

while i<=res,
    EM_($+1) = g^3*(2/9*theta1p(i)*theta2p(i)*cos(theta1(i)-theta2(i))+0.5*theta1p(i)^2-21/18*cos(theta1(i))+4/81*theta2p(i)^2-2/9*cos(theta2(i)))-EM_0;
    i=i+1; end

plot2d(tempo, EM_, 5)
title ("Energia Mecânica em Função do Tempo - M1 & M2")
xlabel("t (s)")
ylabel("EM (J.m/kg)")
legends(leg="M2", 5,opt=3)

show_window(2)

plot2d(tempo, EM_, 5)
title ("Energia Mecânica em Função do Tempo - M2")
xlabel("t (s)")
ylabel("EM (J.m/kg)")

show_window(3)

plot2d(tempo, EM_-EM, 5)
title ("Energia Mecânica em Função do Tempo - M1 & M2")
xlabel("t (s)")
ylabel("EM (J.m/kg)")

