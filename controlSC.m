clc;clear;close all;
tf = 10; % final time
% initial states
rollInt=-0.5; % degree
pitchInt=0.7; % degree
yawInt=1; % degree
rollDotInt=0.1; % degree/sec
pitchDotInt=0.5; % degree/sec
yawDotInt=-0.2; % degree/sec
% Integration
[tn,xn] = ode45(@NonLinSC,[0,10],[-0.5 0.7 1 0.1 0.5 -0.2]);
[tl,xl] = ode45(@LinSC,[0,10],[-0.5 0.7 1 0.1 0.5 -0.2]);
% Plotting
figure(1)
subplot(2,1,1)
plot(tn,xn(:,1:3)); % angles
title('Angles for Nonlinear System','Interpreter','LaTex');
xlabel('Time (sec)','Interpreter','LaTex');
ylabel('Angle (deg)','Interpreter','LaTex');
legend('Roll, $\Phi$','Pitch, $\theta$','Yaw, $\Psi$','Interpreter','LaTex');
grid on;
subplot(2,1,2)
plot(tn,xn(:,4:6)); % angular rates
title('Angular Rates for Nonlinear System','Interpreter','LaTex');
xlabel('Time (sec)','Interpreter','LaTex');
ylabel('Angular rate (deg/sec)','Interpreter','LaTex');
legend('Roll rate, $\dot{\Phi}$','Pitch rate, $\dot{\theta}$','Yaw rate, $\dot{\Psi}$','Interpreter','LaTex');
grid on;
figure(2)
subplot(2,1,1)
plot(tl,xl(:,1:3)); % angles
title('Angles for Linearised System','Interpreter','LaTex');
xlabel('Time (sec)','Interpreter','LaTex');
ylabel('Angle (deg)','Interpreter','LaTex');
legend('Roll, $\Phi$','Pitch, $\theta$','Yaw, $\Psi$','Interpreter','LaTex');
grid on;
subplot(2,1,2)
plot(tl,xl(:,4:6)); % angular rates
title('Angular Rates for Linearised System','Interpreter','LaTex');
xlabel('Time (sec)','Interpreter','LaTex');
ylabel('Angular rate (deg/sec)','Interpreter','LaTex');
legend('Roll rate, $\dot{\Phi}$','Pitch rate, $\dot{\theta}$','Yaw rate, $\dot{\Psi}$','Interpreter','LaTex');
grid on;
figure(3)
subplot(2,1,1)
plot(xn(:,1),xn(:,4)); % phase plane
hold on;
plot(xn(:,2),xn(:,5)); % phase plane
plot(xn(:,3),xn(:,6)); % phase plane
title('Phase Plane for Nonlinear System','Interpreter','LaTex');
xlabel('Angle (deg)','Interpreter','LaTex');
ylabel('Angular rate (deg/sec)','Interpreter','LaTex');
legend('Roll, $\Phi$','Pitch, $\theta$','Yaw, $\Psi$','Interpreter','LaTex');
grid on;
hold off;
subplot(2,1,2)
plot(xl(:,1),xl(:,4));
hold on;
plot(xl(:,2),xl(:,5));
plot(xl(:,3),xl(:,6));
grid on;
title('Phase Plane for Linearised System','Interpreter','LaTex');
xlabel('Angle (deg)','Interpreter','LaTex');
ylabel('Angular rate (deg/sec)','Interpreter','LaTex');
legend('Roll, $\Phi$','Pitch, $\theta$','Yaw, $\Psi$','Interpreter','LaTex');
hold off;
% LQR Controller Code
function K = ControlDesign()
mu = 398600.4418e9;
r = 6378.14e3 + 500e3;
w0 = sqrt(mu/r^3);
J1 = 1700;
J2 = 2000;
J3 = 1400;
a1 = 1/J1*(3*w0^2*(J3-J2));
a2 = 1/J1*(w0*(J2-J3));
a3 = 1/J2*(2*w0^2*(J3-J1));
a4 = 1/J3*(w0*(J1-J2));
A = [0 0 -w0 1 0 0;
    0 0 0 0 1 0;
    w0 0 0 0 0 1;
    a1 0 0 0 0 a2;
    0 a3 0 0 0 0;
    0 0 0 a4 0 0];
B = [zeros(3);
    [1/J1 0 0;0 1/J2 0;0 0 1/J3]];
K = place(A,B,[-1 -2 -3 -4 -5 -6]);
end
% Nonlinear System Dynamics
function deriv = NonLinSC(t,x)
J1 = 1700;
J2 = 2000;
J3 = 1400;
Omega = [x(4);x(5);x(6)];
K = ControlDesign();
u = -K*x;
phidot = Omega(1);
thetadot = Omega(2);
psidot = Omega(3);
Omega1dot = 1/J1*((J2-J3)*Omega(2)*Omega(3)+u(1));
Omega2dot = 1/J2*((J3-J1)*Omega(1)*Omega(3)+u(2));
Omega3dot = 1/J3*((J1-J2)*Omega(2)*Omega(1)+u(3));
deriv=[phidot;thetadot;psidot;Omega1dot;Omega2dot;Omega3dot;];
end
% Linear System Dynamics
function deriv = LinSC(t,x)
mu = 398600.4418e9;
r = 6378.14e3 + 500e3;
w0 = sqrt(mu/r^3);
J1 = 1700;
J2 = 2000;
J3 = 1400;
Omega = [x(4);x(5);x(6)];
K = ControlDesign();
u = -K*x;
phidot = Omega(1)-w0*x(3);
thetadot = Omega(2);
psidot = Omega(3)+w0*x(1);
Omega1dot = 1/J1*(w0*(J2-J3)*Omega(3)*3*w0^2*(J3-J2)*x(1)+u(1));
Omega2dot = 1/J2*(3*w0^2*(J3-J1)*x(2)+u(2));
Omega3dot = 1/J3*(w0*(J1-J2)*Omega(1)+u(3));
deriv=[phidot;thetadot;psidot;Omega1dot;Omega2dot;Omega3dot;];
end
