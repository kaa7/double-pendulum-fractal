clc; close all; clear;
% Variabiles
g = 9.8067;
m1 = 1;
m2 = 2;
l1 = 1;
l2 = 2;
miu = (m1 + m2) / m2;
r = l1 / l2;
w1s = g / l1;
w2s = g / l2;

% Initial conditions
N1 = 10;
N2 = 10;
theta1v = linspace(-pi, pi, N1);
theta2v = linspace(-pi, pi, N2);
culoare = zeros(N1,N2);
for j1 = 1 : N1
    for j2 = 1 : N2
        % Time and preinitialisation
        theta10 = theta1v(j1);
        theta20 = theta2v(j2);
        t0 = 0;
        tf = 5;
        N = 1000;
        t = linspace(t0, tf, N);
        dt = t(2) - t(1);
        dts = dt * dt;
        theta1 = zeros(1,N);
        theta2 = zeros(1,N);
        theta1(1) = theta10;
        theta2(1) = theta20;
        theta1(2) = theta1(1);
        theta2(2) = theta2(1);
        A = [miu * r, 0; 0, 1/r];
        b = zeros(2,1);
        eps = zeros(2,1);
        % Simulate
        for i = 2 : N - 1
            delta = theta2(i) - theta1(i);
            A(1,2) = cos(delta);
            A(2,1) = A(1,2);
            vu1 = (theta1(i) - theta1(i - 1)) / dt;
            vu2 = (theta2(i) - theta2(i - 1)) / dt;
            b = [vu2^2 * sin(delta) - miu * w2s  * sin(theta1(i)); -vu1^2 * sin(delta) - w1s  * sin(theta2(i))];
            eps = A \ b; % eps = A^-1 * b
            theta1(i + 1) = 2 * theta1(i) - theta1(i - 1) + dts * eps(1);
            theta2(i + 1) = 2 * theta2(i) - theta2(i - 1) + dts * eps(2);
        endfor

        t1f = theta1(N);
        t2f = theta2(N);
        culoare(j1,j2) = sin(t1f)*sin(t2f);
    endfor
endfor
surf(theta1v, theta2v, culoare);
grid off
axis off
box off
%shading interp;
view([0,90]);
print -dpdf fractal.pdf;
