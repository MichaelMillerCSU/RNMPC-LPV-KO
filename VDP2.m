clc
clear
close all
r = randi(5000);
rng(3506) 
f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
n = 2;
m = 1; % number of control inputs


A = {};
B = {};
C = {};
Ylift = {};
Xlift = {};
U_sub = {};
X_sub = {};
W = {};
V = {};
VVt = {};
WVt = {};
M = {};



%% ************************** Discretization ******************************

deltaT = 0.05;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Basis functions *****************************
basisFunction = 'rbf';
Nrbf = 2;
cent = rand(2,Nrbf)*4 - 2;
% cent = [[0.7423; 0.9025] [-0.8372; .2427]]
rbf_type = 'gauss'; 
% Lifting mapping - RBFs + the state itself
% liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] - [zeros(2, 1);rbf(zeros(2, 1),cent,rbf_type)]);
liftFun = @(xx)( [xx]);

Nsim = 100;
Ntraj = 100;
Nd = 2;

Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 2*rand([Nsim Ntraj]) - 1;
Xcurrent = (rand(n,Ntraj)*4 - 2);

N = Nsim * Ntraj / Nd;

[AA, BB, CC, YYlift, XXlift, XX, UU] = VDP_System_Generate(liftFun, Ubig, Xcurrent, Nsim, Ntraj);

for i = 1 : Nd
    Ylift{i} = YYlift(:, N * (i - 1) + (1 : N));
    Xlift{i} = XXlift(:, N * (i - 1) + (1 : N));
    X_sub{i} = XX(:, N * (i - 1) + (1 : N));
    U_sub{i} = UU(:, N * (i - 1) + (1 : N));
    W{i} = [Ylift{i} ; X_sub{i}];
    V{i} = [Xlift{i} ; U_sub{i}];
    VVt{i} = V{i}*V{i}';
    WVt{i} = W{i}*V{i}';
    M{i} = WVt{i} * pinv(VVt{i});
    
    A{i} = M{i}(1:Nlift,1:Nlift);
    B{i} = M{i}(1:Nlift,Nlift+1:end);
    C{i} = M{i}(Nlift+1:end,1:Nlift);
end

Nd_temp = 2 * Nd;

for i = Nd + 1 : Nd_temp
    A{i} = -A{i - Nd};
    B{i} = -B{i - Nd};
    C{i} = C{i - Nd};
end

A{Nd_temp + 1} = diag(0.001 * ones(Nlift, 1));
B{Nd_temp + 1} = 0.001 * ones(Nlift, m);
C{Nd_temp + 1} = C{Nd_temp};

Nd = Nd_temp + 1;



%% Not polytopic, just one model
notpolytopic = 0;

if notpolytopic == 1
    for i = 1 : Nd
        A{i} = AA;
        B{i} = BB;
    end
    epsilon = [0.0; 0.0; 0.0; 0.0; 0.0; zeros(Nd - 5, 1)];
end



%% LQR
Q1 = diag([10 * ones(n, 1); 0.1 * ones(Nlift - n, 1)]);
% Q1 = 100 * eye(Nlift);

R = 0.001;
F_Receding = -dlqr(AA, BB, Q1, R);


%% Verification 
N = 1200;
% x0_init = [-0.06; 0.83];
% x0_init = [-1.3; 1.8];
% x0_init = [0.3; 0.83];
% epsilon = [0.00; 0.0; zeros(Nd - 2, 1)];
if notpolytopic == 0
    epsilon = [0.1 ; 0.0; 0.1; 0.0; 0.0; 0.0 * ones(Nd - 5, 1)];
%     epsilon = [0.0 ; 0.0; 0.0; 0.0; 0.0; 0.0 * ones(Nd - 5, 1)];
end

flag = 1;
if flag == 1
    for j = 1 : Nd 
        A{j} = (1 + epsilon(j)) * A{j};
        B{j} = (1 + epsilon(j)) * B{j};
    end
end

u = 2 * rand(1, N) - 1;
for k = 1 : 1
%     x0_init = 2 * rand(2, 1) - 1;
    x0_init = [-1.3; 1.8];
%     epsilon = [0.01; 0.0; 0.01; -1.01; -1.02; -1.02; 0.02; -1.02; 0.01; 0.01; 0.02; -1.02; 0.02; -1.02; 0.02; 0.02];
%     epsilon = 0.5 * epsilon;
    Cy = [1 0]; 
    
    NSpan = 1 : N;
    x0 = x0_init;
    % u = 2 * rand(1, N) - 1;
    % u = prbs(15, N);
    X = {};
    x = {};
    x_temp = {};
    for i = 1 : Nd
        X{i} = [];
        x{i} = liftFun(x0);
    end
    
    x_rec = liftFun(x0);
    Coefficient_Set = [];
    
    X_Real = [x0];
    X_Recon = [x_rec];
    
    P_i = [];
    
    
    
    for i = 1 : N
        for j = 1 : Nd
            x_candidate{j} = C{j} * (   A{j} * x_rec + B{j}  * u(i)   );
        end
    
        M = cell2mat(x_candidate);
        x0 = f_ud(0, x0, u(i));
        p_i = lsqlin(M, x0,[ [-eye(Nd, Nd) ; ];], [zeros(Nd, 1);],  ones(1, Nd), 1)
%         p_i = lsqlin(M, x0, [-eye(Nd, Nd) ; ones(1, Nd)], [zeros(Nd, 1); 1]);
        Coefficient_Set = [Coefficient_Set p_i];
    
        A_Nd = cell2mat(A); 
        B_Nd = cell2mat(B);
    
        A_Rec = A_Nd * kron(p_i, eye(Nlift));
        B_Rec = B_Nd * p_i;
    
        x_rec =  A_Rec * x_rec + B_Rec * u(i);
        
        X_Real = [X_Real x0];
        X_Recon = [X_Recon x_rec];
    end
    
    figure
    plot(X_Real(1, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on 
    plot(X_Recon(1, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
    legend('Original', 'LPV')
    
    figure
    plot(X_Real(2, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on 
    plot(X_Recon(2, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
    legend('Original', 'LPV')
    pause(0.1)
end

%% Plot polytopic coefficient 
figure
for i = 1 : Nd
    plot(Coefficient_Set(i, :), 'LineStyle','--')
    hold on
end


%%
Control_Bound = 4;

Loops = 1;
X_Loop = [];
U_Loop = [];
for L = 1 : Loops

% A1 = A(0.1);
% A2 = A(10);
nx = size(A{1}, 1);
nu = size(B{1}, 2);
init = [1 ; -0.5];
x0 = liftFun(init);
% Q1 = 1 * C1'*C1;

% 
% F_Receding = F_Receding_Static;
X_Receding = [init];
F_Set = [];
x_next = init;
U_Receding = [];
time = 200;
J_Set = [0];
J  = 0;
for i = 1 : time
    gamma = sdpvar(1, 1);
    X = sdpvar(nu, nu);
    Q = sdpvar(nx, nx);
    Y = sdpvar(nu, nx);
    LMI1 = [1 liftFun(x_next)';
            liftFun(x_next) Q];

    for j = 1 : Nd
        LMI2{j} = [Q                       (A{j} * Q + B{j} * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
                   A{j} * Q + B{j} * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
                   sqrt(Q1)*Q              zeros(nx, nx)        (gamma) * eye(nx, nx) zeros(nx, nu);
                   sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       (gamma) * eye(nu, nu)];
    end

    LMI0 = [X  Y;
            Y' Q];

    ZEROS = 0;
    Constraints = [LMI0 >= 0, LMI1 >= ZEROS, Q >= ZEROS];
    for j = 1 : Nd
        Constraints = [Constraints LMI2{j} >= ZEROS];
    end

    % Input constraints  |uk| <= 2

    for j = 1 : nu
        Constraints = [Constraints X(j, j) <= Control_Bound^2 ] ;
    end

    Objective = gamma;

    sol  = solvesdp(Constraints, Objective)
    if sol.problem ~= 0
        error("Infeasible!")
    end
    Y = double(Y)
    Q = double(Q)
    
    F_Receding = Y / Q
        
    J = J + x_next' * x_next + [F_Receding * liftFun(x_next)]'*[F_Receding * liftFun(x_next)]
    J_Set = [J_Set J];
    x_next = f_ud(0,x_next,F_Receding * liftFun(x_next))
%     x_next = f_ud(0,x_next,0)
    X_Receding = [X_Receding x_next];
    U_Receding = [U_Receding   F_Receding * liftFun(x_next)];
    F_Set  = [F_Set; F_Receding];
end
    U_Loop = [U_Loop U_Receding];
    X_Loop = [X_Loop X_Receding];
end


Tspan = 0:0.05: time * 0.05;
NTime = size(Tspan, 2);
figure
subplot 211
for L = 0 : Loops - 1
    plot(Tspan,X_Loop(1, L * (time + 1) + 1 : (L + 1)* (time + 1)), 'LineWidth', 3, 'Color', 'b')
    hold on 
    plot(Tspan,X_Loop(2, L * (time + 1) + 1 : (L + 1)* (time + 1)), 'LineWidth', 3, 'Color', 'r')
end
% axis([0 10,-0.2 1]);
xlabel('$Time(sec)$','interpreter','latex');
ylabel('$State$','interpreter','latex');
title('$State\,\,trajectory$','interpreter','latex');
legend('$x_1$', '$x_2$', 'Interpreter', 'latex')


subplot 212
for L = 0 : Loops - 1
    plot(Tspan(1: end - 1),U_Loop(1, L * (time + 1) + 1 : (L + 1)* (time)), 'LineWidth', 3)
    hold on 
end
% plot(Tspan,Control_Bound * ones(1, NTime),'LineWidth', 3, 'LineStyle','--');
% plot(Tspan,-Control_Bound* ones(1, NTime),'LineWidth', 3, 'LineStyle','--');

% axis([0 10,-2 0.5]);
xlabel('time(sec)');
ylabel('$u$ ','interpreter','latex');
title('Control signal $u$','interpreter','latex');




figure
plot(Tspan(1 : end), J_Set,'LineWidth', 3, 'LineStyle','-')
xlabel('time(sec)');
ylabel('$J$ ','interpreter','latex');





