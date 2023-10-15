// Sistema de EDOs Linear

function dxdt = f(t,x)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    dxdt = [x3; x4; (-7/3)*x1+x2; (21/4)*x1-(9/2)*x2];
endfunction

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

// Gráfico da Trajetória Linear (M1)

show_window(0)

g = 9.81 // Valor de convenção em m/s^2
theta1 = S_x(1,:);
theta2 = S_x(2,:);

x1data = g*sin(theta1);
y1data = -g*cos(theta1);

x2data = x1data+2/3*g*sin(theta2);
y2data = y1data-2/3*g*cos(theta2);

plot(x2data,y2data,"r")
legends(leg="Modelagem Linear",5,opt=3)
a=get("current_axes")

// Solução do Sistema Não-Linear (M1)

SN_x = ode("rk",x0,t0,tempo,h);

// Gráfico da Trajetória Não-Linear (M1)

theta1N = SN_x(1,:);
theta2N = SN_x(2,:);

x1dataN = g*sin(theta1N);
y1dataN = -g*cos(theta1N);

x2dataN = x1dataN+2/3*g*sin(theta2N);
y2dataN = y1dataN-2/3*g*cos(theta2N);

plot(x2dataN,y2dataN,"b")
xlabel("x (m)")
ylabel("y (m)")
title ("Trajetória do Extremo Livre - C1, M1")
legends(leg="Modelagem Não-Linear",2,opt=4)
a=get("current_axes")

// Solução do Sistema Linear (M2)

S_x = ode("adams",x0,t0,tempo,f);

// Gráfico da Trajetória Linear (M2)

show_window(1)

theta1 = S_x(1,:);
theta2 = S_x(2,:);

x1data = g*sin(theta1);
y1data = -g*cos(theta1);

x2data = x1data+2/3*g*sin(theta2);
y2data = y1data-2/3*g*cos(theta2);

plot(x2data,y2data,"r")
legends(leg="Modelagem Linear", 5,opt=3)

// Solução do Sistema Não-Linear (M2)

SN_x = ode("adams",x0,t0,tempo,h);

// Gráfico da Trajetória Não-Linear (M2)

theta1N = SN_x(1,:);
theta2N = SN_x(2,:);

x1dataN = g*sin(theta1N);
y1dataN = -g*cos(theta1N);

x2dataN = x1dataN+2/3*g*sin(theta2N);
y2dataN = y1dataN-2/3*g*cos(theta2N);

plot(x2dataN,y2dataN,"b")
xlabel("x (m)")
ylabel("y (m)")
title ("Trajetória do Extremo Livre - C1, M2")
legends(leg="Modelagem Não-Linear",2,opt=4)
