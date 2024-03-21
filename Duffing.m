clc
clear
close all
r = randi(5000);
rng(3506) 
f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
n = 2;
m = 1; % number of control inputs


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
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] - [zeros(2, 1);rbf(zeros(2, 1),cent,rbf_type)]);

Nsim = 200;
Ntraj = 100;

Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 2*rand([Nsim Ntraj]) - 1;
Xcurrent = (rand(n,Ntraj)*4 - 2);

N = Nsim * Ntraj / 2;

[A, B, C, Ylift, Xlift, X, U] = Duffing_System_Generate(liftFun, Ubig, Xcurrent, Nsim, Ntraj);

Ylift1 = Ylift(:, (1 : N));
Xlift1 = Xlift(:, (1 : N));
X1 = X(:, (1 : N));
U1 = U(:, (1 : N));
W1 = [Ylift1 ; X1];
V1 = [Xlift1; U1];
V1Vt = V1*V1';
W1Vt = W1*V1';
M1 = W1Vt * pinv(V1Vt);

A1 = M1(1:Nlift,1:Nlift);
B1 = M1(1:Nlift,Nlift+1:end);
C1 = M1(Nlift+1:end,1:Nlift);


Ylift2 = Ylift(:, N + (1 : N));
Xlift2 = Xlift(:, N + (1 : N));
X2 = X(:, N + (1 : N));
U2 = U(:, N + (1 : N));
W2 = [Ylift2 ; X2];
V2 = [Xlift2; U2];
V2Vt = V2*V2';
W2Vt = W2*V2';
M2 = W2Vt * pinv(V2Vt);

A2 = M2(1:Nlift,1:Nlift);
B2 = M2(1:Nlift,Nlift+1:end);
C2 = M2(Nlift+1:end,1:Nlift);


X11 = [];
X22 = [];
for i = 1 : Ntraj : N
    X11 = [X11 X1(:, i)];
    X22 = [X22 X2(:, i)];
end

plot([X11(1, 1:N / Ntraj) X22(1, 1)], [X11(2,1:N / Ntraj)  X22(2,1)], 'LineWidth', 2.0)
hold on 
plot(X22(1, 1:N / Ntraj), X22(2,1:N / Ntraj), 'LineWidth', 2.0)
% norm(A3 - A4)

%% LQR
Q1 = diag([[10 10 0.1 0.1]]);
% Q1 = 100 * eye(Nlift);

R = 0.001;
F_Receding = -dlqr(A, B, Q1, R);


%% Verification 
N = 1200;
% x0_init = [-0.06; 0.83];
% x0_init = [-1.3; 1.8];
x0_init = [0.3; 0.83];

% x0_init = 2 * rand(2, 1) - 1;
u = 2 * rand(1, N) - 1;
epsilon = 0.01;
Cy = [1 0]; 
for j = 1 : 2
    NSpan = 1 : N;
    x0 = x0_init;
    % u = 2 * rand(1, N) - 1;
    % u = prbs(15, N);
    X1 = [];
    X2 = [];
    X3 = [];
    X4 = [];
    X = [];
    X_Recon = [];
    
    x1 = liftFun(x0);
    x2 = liftFun(x0);
    x3 = liftFun(x0);
    x4 = liftFun(x0);
    x_rec = liftFun(x0);
    a_b_Set = [];
    if j == 2 
        Delta_p = [];
    end
    x_next = x0;
    
    for i = 1 : N
        X1 = [X1 x1];
        X2 = [X2 x2];
        X = [X x0];
        X_Recon = [X_Recon x_rec];
    
        x_1 = C1 * A1 * x_rec + C1 *  B1 * u(i);
        x_2 = C2 * A2 * x_rec + C2 *  B2 * u(i);
        xlift = liftFun(x0);
    
        x0 = f_ud(0, x0, u(i));
        U0 = u(i);
        ylift = liftFun(x0);
    
        if j == 1
            a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1)
        end
        if j == 2
            a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ;  ones(1, 2)];], [zeros(2, 1); 1 ])

        end
        a_b_Set = [a_b_Set a_b]; 
        
        x1 = A1 * x1 + B1 * u(i);
        x2 = A2 * x2 + B2 * u(i);
        A = a_b(1) * A1 + a_b(2) * A2;
        B = a_b(1) * B1 + a_b(2) * B2;
        if j == 2
            delta_p = 1 - sum(a_b_Set(:, i));
            Delta_p = [Delta_p delta_p];
            A = a_b(1) * A1 + a_b(2) * A2 ;
            B = a_b(1) * B1 + a_b(2) * B2 ;
        end
        
        x_rec =  A * x_rec + B * u(i);
    
    
    end
    
    

if j == 1
    figure
    subplot 311
%     plot(NSpan, X1(1, :))
%     hold on
%     plot(NSpan, X2(1, :))
%     hold on
    plot(NSpan * deltaT,X(1, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on
    plot(NSpan * deltaT, X_Recon(1, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
    legend('Original System','LPV System', 'Interpreter', 'latex')
    ylabel('State $x_1$',  'Interpreter', 'latex')

    subplot 312
%     plot(NSpan, X1(2, :))
%     hold on
%     plot(NSpan, X2(2, :))
%     hold on
    plot(NSpan * deltaT,X(2, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on
    plot(NSpan * deltaT, X_Recon(2, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
    legend('Original System','LPV System', 'Interpreter', 'latex')
    ylabel('State $x_2$',  'Interpreter', 'latex')
    
    subplot 313
    plot(NSpan * deltaT, a_b_Set(1, :))
    hold on
    plot(NSpan * deltaT, a_b_Set(2, :))
    legend('$p_1$', '$p_2$', 'Interpreter', 'latex')
    xlabel('$Time(sec)$',  'Interpreter', 'latex')
    ylabel('Coefficient $p_i$',  'Interpreter', 'latex')

end
    if j == 2
        figure
        subplot 411
        plot(NSpan * deltaT,X(1, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
        hold on
        plot(NSpan * deltaT, X_Recon(1, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
        legend('Original System','LPV System', 'Interpreter', 'latex')
        ylabel('State $x_1$',  'Interpreter', 'latex')
    
        subplot 412
        plot(NSpan * deltaT,X(2, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
        hold on
        plot(NSpan * deltaT, X_Recon(2, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
        legend('Original System','LPV System', 'Interpreter', 'latex')
        ylabel('State $x_2$',  'Interpreter', 'latex')
        
        subplot 413
        plot(NSpan * deltaT, a_b_Set(1, :))
        hold on
        plot(NSpan * deltaT, a_b_Set(2, :))
        legend('$p_1$', '$p_2$', 'Interpreter', 'latex')
        ylabel('Coefficient $p_i$',  'Interpreter', 'latex')

        subplot 414
        plot(NSpan * deltaT, Delta_p, 'LineWidth', 2, 'LineStyle', '-','Color','g')
        legend('$\Delta_p$', 'Interpreter', 'latex')
        ylabel('Coefficient $p_0/\Delta_p$',  'Interpreter', 'latex')

        xlabel('$Time(sec)$',  'Interpreter', 'latex')

    end

    pause
    close all

    if j == 1
    A1 = (1 + epsilon) * A1;
    B1 = (1 + epsilon) * B1;
    end
end


%% Figure 5 Closed-loop responses for the time-varying system with input constraint; 

Control_Bound = 30;

Loops = 1;
X_Loop = [];
U_Loop = [];
for L = 1 : Loops

% A1 = A(0.1);
% A2 = A(10);
nx = size(A1, 1);
nu = size(B1, 2);
init = [-1 ; 1];
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

    LMI2 = [Q                       (A1 * Q + B1 * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
            A1 * Q + B1 * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
            sqrt(Q1)*Q              zeros(nx, nx)        (gamma) * eye(nx, nx) zeros(nx, nu);
            sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       (gamma) * eye(nu, nu)];
    
    LMI3 = [Q                       (A2 * Q + B2 * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
            A2 * Q + B2 * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
            sqrt(Q1)*Q              zeros(nx, nx)        (gamma) * eye(nx, nx) zeros(nx, nu);
            sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       (gamma ) * eye(nu, nu)];


    LMI0 = [X  Y;
            Y' Q];

    ZEROS = 0;
    Constraints = [LMI0 >= 0, LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, Q >= ZEROS];

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
        
%     f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );
    %     x_next = A_Receding * x_next + B * F_Receding * x_next;
    Xlift1 = [Xlift1 liftFun(x_next)];
    Xlift2 = [Xlift2 liftFun(x_next)];
%     Xlift3 = [Xlift3 liftFun(x_next)];
%     Xlift4 = [Xlift4 liftFun(x_next)];
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





