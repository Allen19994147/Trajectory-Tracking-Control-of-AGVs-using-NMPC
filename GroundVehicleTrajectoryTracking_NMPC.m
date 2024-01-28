%%% Ground Vehicle Trajectory Tracking Control using NMPC
%%% Allen Lee
clear all
clc



%% Model Definition
global timestep
global N n m x xr ur
N = 20; % prediction horizen
n = 3; % # state
m = 2; % # input
x = zeros(n,1); % vehicle pose
u = zeros(m,1); % velocity of left&right wheels

timestep = 0.01;
Simtime = 10; % 20 seconds
timeline = 0:timestep:Simtime-timestep;
RunSample = Simtime/timestep;
X = zeros(n,RunSample+N+1); U = zeros(m,RunSample+N);
Xr = zeros(n,size(X,2)+N+1); Ur = zeros(m,size(U,2)+N); % Equilibrium points

% Trajectory Definition 
Xr(1,N+2:RunSample+2+N) = (0:RunSample).*timestep*0.5;
Xr(2,N+2:RunSample+2+N) = 0.5.*sin(0.5.*(0:RunSample).*timestep);
Xr(1:2,RunSample+3+N:end) = ones(1,size(Xr,2)-RunSample-N-2).*Xr(1:2,RunSample+N+1);

% figure
% scatter(Xr(1,:) ,Xr(2,:))
% figure
% hold on
% plot(Xr(1,:))
% plot(Xr(2,:))
% hold off
%% Trajectory Tracking using NMPC
xr = zeros(n,N);
ur = zeros(m,N);

Xk = zeros(n*(N+1),1); % all future x from step k
Uk = zeros(m*N,1); % all future u from step k
Xk(1:n,1) = x;
Uk(1:m,1) = u;

zek = [Xk;Uk];
% Difference from the references
zek(1:n,1) = x-xr(:,1); 
zek((N+1)*n+1:(N+1)*n+m) = u-ur(:,1);

% Constraint Definition
xmin = [0;-0.501;-pi];
xmax = [10.01;0.501;pi];
umin = [-80;-80];
umax = [80;80];

Fx = [eye(n);-eye(n)];
Fu = [eye(m);-eye(m)];

R = eye(m);
Qx = diag([1 1 0]);
QX = zeros(n*(N+1),n*(N+1));
RU = zeros(m*N,m*N);
FX = zeros(2*n*(N+1),n*(N+1));
FU = zeros(2*m*N,m*N);
Gxe = zeros(2*n*(N+1),1);
Gue = zeros(2*m*N,1);
for i = 1:N
    QX((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Qx;
    RU((i-1)*m+1:i*m,(i-1)*m+1:i*m) = R;

    FX((i-1)*2*n+1:i*2*n,(i-1)*n+1:i*n) = Fx;
    FU((i-1)*2*m+1:i*2*m,(i-1)*m+1:i*m) = Fu;

end

QX(N*n+1:end,N*n+1:end) = Qx;
FX(N*2*n+1:end,N*n+1:end) = Fx;
H = blkdiag(QX,RU.*0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with NMPC
ObjFunc = @(ze) (ze'*H*ze);
F = blkdiag(FX,FU);
Feq = []; Geq = []; lb = []; ub = []; % Linear constraints
nonlcon = @NonlinearConstraints;
options = optimoptions('fmincon','Display','none','Algorithm','sqp');

for i = 1:RunSample+N
    x = X(:,i); % for nonlon
    xr(:,1:end-1) = xr(:,2:end);
    xr(:,end) = Xr(:,i);
    ur(:,1:end-1) = ur(:,2:end);
    ur(:,end) = Ur(:,i);
    
    Gxe = zeros(2*n*(N+1),1);
    Gue = zeros(2*m*N,1);
    for j = 1:N
        Gxe((j-1)*2*n+1:j*2*n,1) = [xmax;-xmin] - Fx*(xr(:,j));
        Gue((j-1)*2*m+1:j*2*m,1) =  [umax;-umin] - Fu*(ur(:,j));
    end
    Gxe(N*2*n+1:end,1) = [xmax;-xmin] - Fx*(xr(:,end));
    G = [Gxe;Gue];

    zek = fmincon(ObjFunc,zek,F,G,Feq,Geq,lb,ub,nonlcon,options);
    U(:,i) = zek(n*(N+1)+1:n*(N+1)+m,1) + ur(:,1);

    X(:,i+1) =  VehicleDynamic(X(:,i),U(:,i));
end

%% Save and Compare Data
[X(1:2,:)',Xr(1:2,1:size(X,2))'];
rmse(X(1:2,:)',Xr(1:2,1:size(X,2))');
save("TrajectoryTracking_Data.mat","X","Xr","U")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results and Visulization
figure (1)
hold on
scatter(X(1,:),X(2,:))
scatter(Xr(1,:),Xr(2,:))
hold off

figure (2)
hold on
plot(U(1,:))
plot(U(2,:))
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq] =NonlinearConstraints(z)
global N n m x xr ur
c=[];
shfu = (N+1)*n;
for i=1:N
    x_dummy = z((i-1)*n+1:i*n,1)+xr(:,i);
    u_dummy = z(shfu+(i-1)*m+1:shfu+i*m)+ur(:,i);
    x_prime = VehicleDynamic(x_dummy,u_dummy);
    c1((i-1)*n+1:i*n,1)=z(i*n+1:(i+1)*n,1)+xr(:,i) - x_prime;

end
ceq = [z(1:n)-x+xr(:,i);c1];
end



function x_prime = VehicleDynamic(x,u)
global timestep

veh_len = 1.31; % vehicle length
v = 0.5*(u(1)+u(2));
w = (u(1)-u(2))/veh_len;

x_change = [v*timestep*cos(x(3)+w*timestep);
            v*timestep*sin(x(3)+w*timestep);
            w*timestep];
x_prime = x + x_change;

end