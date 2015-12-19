function [] = arenstorf()
clc;

tmax = 17.3;
tspan = [0. tmax];
y0 = [0.994 0. 0. -2.00158510637908 0. 0.];
mu = (7.349/597.4);
options = odeset('RelTol', 1e-5);

[T,Y] = ode45(@(t,y)dgl(t,y,mu),tspan,y0,options);
plot3(cos(2*pi*T/tmax),sin(2*pi*T/tmax),(0.*T));
hold on
plot3((0.*T),(0.*T),(0.*T),'*');
plot3(Y(:,1),Y(:,3),Y(:,5));
hold off
axis equal
end

function [dy] = dgl(t,y,mu)
dy = [0. 0. 0. 0. 0. 0.]';
% size(y)
dy(1) = y(2);
dy(2) = y(1) + 2 * y(4) - ((1 - mu) * (y(1) + mu))/(((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + ...
    y(5)*y(5)) * sqrt((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + y(5)*y(5))) - ...
    (mu * (y(1) - 1 + mu))/(((y(1) + mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)) * ...
    sqrt((y(1) + mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)));
dy(3) = y(4);
dy(4) = y(3) - 2 * y(2) - ((1 - mu) * y(3))/(((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + ...
    y(5)*y(5)) * sqrt((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + y(5)*y(5))) - (mu * ...
    y(3))/(((y(1) + mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)) * sqrt((y(1) + ...
    mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)));
dy(5) = y(6);
dy(6) = -((1 - mu) * y(5))/(((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + y(5)*y(5)) * ...
    sqrt((y(1) + mu)*(y(1) + mu) + y(3)*y(3) + y(5)*y(5))) - (mu * ...
    y(5))/(((y(1) + mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)) * ...
    sqrt((y(1) + mu - 1)*(y(1) + mu - 1) + y(3)*y(3) + y(5)*y(5)));
end