clear all
clc
close all

delete *.fig
delete *.pdf
delete *.mat

global n vm L
global K
global H
global tau amc_max
global Td tf_l eta_l_1 eta_l_2
global tf t1 eta_f_1 eta_f_2

%%%%%%% Follower graph %%%%%%%%%%%%%%%%%%%%%%
G = graph([1,2,3,4],[2,3,4,1]);
L = full(laplacian(G));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = numnodes(G)+1;

%%%%%%% Leader incidence matrix %%%%%%%%%%%%%%
B = zeros(size(L));
B(1,1) = 1;

H = L + B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eig_H = eig(H);
eta_m = 1.0/(eig_H(1));

%%%% Parameters %%%%%
Td = 30;

eta_l_1 = 7;
eta_l_2 = 1.5;

tf_l = 5;

eta_f_1 = 15;
eta_f_2 = 4;

tf = 25;
t1 = tf - 1;
tau = 0.5*ones(n,1);

amc_max = 300;

N = 3;
K = 2*N-1;

%%%% Parameters %%%%%


%%% Initial Conditions %%%%
r0 = 5000*ones(n,1);
vm = 200*ones(n,1);     % vm is constant, may be different for different missiles
theta0 = [45;0;-45;60;120];
theta0 = deg2rad(theta0);

gamma0 = [90;60;-30;30;70];
gamma0 = deg2rad(gamma0);

x0 = -r0.*cos(theta0);
y0 = -r0.*sin(theta0);

vx0 = vm.*cos(gamma0);
vy0 = vm.*sin(gamma0);

am0 = zeros(n,1);

tspan = 0:0.01:Td;
s0 = [x0;y0;gamma0;am0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derv_func(0,s0);
% % output(0,s0);
% close all

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.1);
[t,s] = ode45(@derv_func, tspan, s0, opts);

x = s(:,1:n);
y = s(:,n+1:2*n);
gamma = s(:,2*n+1:3*n);
am = s(:,3*n+1:4*n);

%%% Plots %%%%%%%%%%%%%%
i_plot = 0;
i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'graph';
plot(G,'LineWidth',2);

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'xy';
plot(x,y,'LineWidth',2)
% xlim([-5000,5000]),ylim([-5000,5000])
hold on
figure(i_plot), scatter(x0,y0,'o','filled','MarkerFaceColor','r')
hold on
quiver(x0,y0,vx0,vy0,0.25,'r')
xlabel('X position (m)'), ylabel('Y position (m)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

tgo = zeros(length(t),n);
tgodot = zeros(length(t),n);
r = zeros(length(t),n);
u = zeros(length(t),n);
am = zeros(length(t),n);
amc = zeros(length(t),n);
error_f = zeros(length(t),n);
thetaM = zeros(length(t),n);
for i=1:length(t)
    [a,b,c,d,e,f,g,hh] = output(t(i),s(i,:)');
    tgo(i,:) = a';
    r(i,:) = b';
    u(i,:) = c';
    am(i,:) = d';
    amc(i,:) = e';
    error_f(i,:) = f';
    tgodot(i,:) = g';
    thetaM(i,:) = hh';
end

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'r';
plot(t,r,'LineWidth',2)
xlabel('Time (s)'), ylabel('Radial distance (m)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'tgo';
plot(t,tgo,'LineWidth',2)
xlabel('Time (s)'), ylabel('Time-to-go (s)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on
hold on
plot([0,Td],[Td,0],'--k','HandleVisibility','off')

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'tgodot';
plot(t,tgodot,'LineWidth',2)
xlabel('Time (s)'), ylabel('Derivative of time-to-go ($/s$)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'u';
plot(t,u,'LineWidth',2)
% ylim([-20,20])
xlabel('Time (s)'), ylabel('Auxiliary control input (/s)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'am';
plot(t,am,'LineWidth',2)
xlabel('Time (s)'), ylabel('Lateral acceleration (m/s^2)')
% ylim([-am_max-10,am_max+10])
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'amc';
plot(t,amc,'LineWidth',2)
xlabel('Time (s)'), ylabel('Lateral acceleration command (m/s^2)')
ylim([-amc_max-10,amc_max+10])
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'error_f';
plot(t,error_f,'k','LineWidth',2)
xlabel('Time (s)'), ylabel('Norm of time-to-go error of followers (s)')
% ylim([-am_max-10,am_max+10])
% legend('1','2','3','4','5')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'thetaM';
plot(t,rad2deg(thetaM),'LineWidth',2)
xlabel('Time (s)'), ylabel('Heading angle error ($^{\circ}$)')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Save figure and data %%%%%%%%%%%%%%%%%%%%%%%%

save output.mat n x y gamma t r tgo u am error_f thetaM

for i=1:i_plot
    savefig(h(i),h(i).Name);
%     exportgraphics(h(i),h(i).Name+".pdf",'ContentType','vector');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = derv_func(t,s)
dxdt = fx(t,s);
end

function [tgo,r,u,am,amc,error_f,tgodot,thetaM] = output(t,s)

[dxdt,tgo,r,u,am,amc,error_f,tgodot,thetaM] = fx(t,s);

end

function [dxdt,tgo,r,u,am,amc,error_f,tgodot,thetaM] = fx(t,s)

global n vm L
global K
global H
global tau amc_max
global Td tf_l eta_l_1 eta_l_2
global tf t1 eta_f_1 eta_f_2

t

x = s(1:n);
y = s(n+1:2*n);
gamma = s(2*n+1:3*n);
am = s(3*n+1:4*n);

vx = vm.*cos(gamma);
vy = vm.*sin(gamma);

r = sqrt(x.^2 + y.^2);
theta = atan2(-y,-x);

iComplete = find(r < 5);

thetaM = gamma - theta;

thetaM = mod(thetaM,2*pi);
for i=1:n
    if (thetaM(i) > pi)
        thetaM(i) = thetaM(i) - 2*pi;
    end
end
% thetaM = atan(tan(thetaM));

vr = -vm.*cos(thetaM);
vtheta = -vm.*sin(thetaM);

tgo = (r./vm).*(1 + ((thetaM.^2)/(2*K)));

F = (vr.*tgo)./r - (vtheta.*thetaM)./(vm*K);
B = (r.*thetaM)./(vm.^2*K);
tgodot = F + B.*am;

error = norm(tgo);

tgo_l = tgo(1);
tgo_f = tgo(2:n);

tgodot_l = tgodot(1);
tgodot_f = tgodot(2:n);

%%% for leader
e_l = tgo_l - (Td-t);
edot_l = tgodot_l + 1;

e_l = e_l/100;
edot_l = edot_l/100;

psi1 = (eta_l_1/(tf_l-t))*(1-exp(-e_l));
dpsi1_dt = (eta_l_1/(tf_l-t)^2)*(1-exp(-e_l));
dpsi1_dx = (eta_l_1/(tf_l-t))*(exp(-e_l));

z2 = edot_l + psi1;
psi2 = (eta_l_2/(tf_l-t))*(1-exp(-z2));
if (t < tf_l-1e-6)
    u_l = -e_l - dpsi1_dt - dpsi1_dx*edot_l - psi2;
else 
    u_l = 0;
end
%%% 

%%% for follower
tgod_f = tgo_l*ones(length(tgo_f),1);
tgod_dot_f = -1*ones(length(tgo_f),1);

e_f = (tgo_f - tgod_f);
error_f = norm(e_f);

edot_f = (tgodot_f - tgod_dot_f);

e_f = e_f/100;
edot_f = edot_f/100;

one_n_f = ones(length(e_f),1);

phi = (eta_f_1/(tf-t))*(one_n_f - exp(-H*e_f));
dphi_dx = (eta_f_1/(tf-t))*diag(exp(-H*e_f))*H;
dphi_dt = (eta_f_1/(tf-t)^2)*(one_n_f - exp(-H*e_f));

z = edot_f + phi;

if (t < tf_l)
    u_f = zeros(size(e_f));
elseif (t >= tf_l && t < t1)
    u_f = -dphi_dx*edot_f - dphi_dt - (eta_f_2/(t1-t))*(one_n_f - exp(-z));
elseif (t >= t1 && t < tf - 1e-6)
    u_f = -dphi_dx*edot_f - dphi_dt;   
else 
    u_f = zeros(size(e_f));
end

u = [u_l*100;u_f*100];

Fdot = (-((am.*vtheta)./(vm.*r)) + ((vtheta.^2 - vr.^2)./(r.^2))).*tgo ...
        + (vr./r).*tgodot - ((vr.*thetaM + vtheta)./(vm.^2*K)).*am ...
        + ((vr.*vtheta.*thetaM + vtheta.^2)./(r.*vm*K));

Bdot = (1./(vm.^3*K)).*(thetaM.*vr.*vm + r.*am - vm.*vtheta);

amc = (tau./B).*(u - Fdot - Bdot.*am + (B.*am)./tau);

for i=1:n
    if (abs(amc(i)) > amc_max)
        amc(i) = amc_max*sign(amc(i));
    end
end

gammadot = am./vm;

amdot = (amc - am)./tau;

vx(iComplete) = 0;
vy(iComplete) = 0;
gammadot(iComplete) = 0;
am(iComplete) = 0;

dxdt = [vx;vy;gammadot;amdot];

end
