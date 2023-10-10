clear all
clc
close all

delete *.fig
delete *.pdf
delete *.mat

global n vm
global K
global H tf eta_l eta_f
global am_max

%%%%%%% Follower graph %%%%%%%%%%%%%%%%%%%%%%
G = graph([1,2,3,4],[2,3,4,1]);
L = full(laplacian(G));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Leader incidence matrix %%%%%%%%%%%%%%
B = zeros(size(L));
B(1,1) = 1;

H = L + B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eig_H = eig(H);
eta_m = 1.0/(eig_H(1));

eta_l = 3;
eta_f = eta_m*2;

tf = 25;

n = numnodes(G)+1;
r0 = 5000*ones(n,1);
vm = 200*ones(n,1);     % vm is constant, may be different for different missiles
theta0 = [45;0;-45;60;120];
theta0 = deg2rad(theta0);

gamma0 = [90;60;-30;30;70];
% gamma0 = [90;91;92;93;94];
gamma0 = deg2rad(gamma0);

N = 3;
K = 2*N-1;
am_max = 30*10;

w = 0.5*ones(n,1);

x0 = -r0.*cos(theta0);
y0 = -r0.*sin(theta0);

vx0 = vm.*cos(gamma0);
vy0 = vm.*sin(gamma0);

tspan = 0:0.01:30;
s0 = [x0;y0;gamma0];

% derv_func(0,s0);
% % output(0,s0);
% close all

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.1);
[t,s] = ode45(@derv_func, tspan, s0, opts);

x = s(:,1:n);
y = s(:,n+1:2*n);
gamma = s(:,2*n+1:3*n);

%%%%%%%%%%  Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
r = zeros(length(t),n);
error_f = zeros(length(t),n);
for i=1:length(t)
    [a,b,c,d,e] = output(t(i),s(i,:)');
    tgo(i,:) = a';
    r(i,:) = b';
    u(i,:) = c';
    am(i,:) = d';
    error_f(i,:) = e';
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
plot([0,30],[30,0],'--k','HandleVisibility','off')

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'u';
plot(t,u,'LineWidth',2)
xlabel('Time (s)'), ylabel('Auxiliary control input')
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'am';
plot(t,am,'LineWidth',2)
xlabel('Time (s)'), ylabel('Lateral acceleration command ($m/s^2$)')
ylim([-100,am_max+10])
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4')
grid on

i_plot = i_plot + 1;
h(i_plot) = figure(i_plot); h(i_plot).Name = 'error_f';
plot(t,error_f,'LineWidth',2,'Color','k')
xlabel('Time (s)'), ylabel('Norm of time-to-go error of followers (s)')
% ylim([-am_max-10,am_max+10])
% legend('1','2','3','4','5')
grid on

% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Save figure and data %%%%%%%%%%%%%%%%%%%%%%%%

save output.mat n x y gamma t r tgo u am error_f

for i=1:i_plot
    savefig(h(i),h(i).Name);
%     exportgraphics(h(i),h(i).Name+".pdf",'ContentType','vector');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = derv_func(t,s)
dxdt = fx(t,s);
end

function [tgo,r,u,am,error_f] = output(t,s)

[dxdt,tgo,r,u,am,error_f] = fx(t,s);

end

function [dxdt,tgo,r,u,am,error_f] = fx(t,s)

global n vm
global K
global H tf eta_l eta_f
global am_max


x = s(1:n);
y = s(n+1:2*n);
gamma = s(2*n+1:3*n);

vx = vm.*cos(gamma);
vy = vm.*sin(gamma);

r = sqrt(x.^2 + y.^2);
theta = atan2(-y,-x);

iComplete = find(r < 1e-3);

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

tgo_l = tgo(1);
tgo_f = tgo(2:n);

reg_factor = 0.1;

Td = 30;
tf_l = 5;

e_l = (tgo_l - (Td-t));
e_l = e_l*reg_factor;
if (t < tf_l-1e-6)
    u_l = -(eta_l/(tf_l-t))*(1-exp(-e_l));
else 
    u_l = 0;
end


tgod = tgo_l*ones(length(tgo_f),1);
e_f = (tgo_f - tgod);

error_f = norm(e_f);

e_f = e_f*reg_factor;

one_n_f = ones(length(e_f),1);
if (t < tf_l)
    u_f = zeros(size(e_f));
elseif (t >= tf_l && t < tf-1e-6)
    u_f = -(eta_f/(tf-t))*(one_n_f - exp(-H*e_f));
else 
    u_f = zeros(size(e_f));
end

u = [u_l/reg_factor;u_f/reg_factor];

am = zeros(n,1);
for i=1:n
    am(i) = (vm(i)*vtheta(i)/r(i)) + ((K*vm(i)^2*abs(vr(i))*tgo(i))/(r(i)^2*thetaM(i))) - ((vm(i)^2*K)/(r(i)*thetaM(i))) + ((vm(i)^2*K)/(r(i)*thetaM(i)))*u(i);
%     am(i) = (vm(i)*vtheta(i)/r(i)) - ((K*vm(i)^2*vr(i)*tgo(i))/(r(i)^2*thetaM(i))) - ((vm(i)^2*K)/(r(i)*thetaM(i))) + ((vm(i)^2*K)/(r(i)*thetaM(i)))*u(i);
    if (abs(am(i)) > am_max)
        am(i) = am_max*sign(am(i));
    end
end

% am = zeros(n,1);

gammadot = am./vm;

vx(iComplete) = 0;
vy(iComplete) = 0;
gammadot(iComplete) = 0;
am(iComplete) = 0;

dxdt = [vx;vy;gammadot];

end
