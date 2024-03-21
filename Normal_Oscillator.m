clc
clear
close all
r = randi(5000);
% rng(2553)
% rng(3000)
rng(2954)
% rng(21665 )

% %% *************************** Dynamics ***********************************
f_u =  @(t,x,u)([ x(2,:) ; x(1, :) - x(2, :).^3 + u] );
% f_u =  @(t,x,u)([ x(2,:) ;x(1, :).^2+ 0.15 .* u(1, :).^3  + 0.1* (1 + x(2,:).^2).*u(1, :) + sin(0.1 .* u(1, :))] );
% f_u =  @(t,x,u)([ x(2,:) ; 4*9.8*sin(x(1, :)) - u .* cos(x(1, :))] );
% f_u =  @(t,x,u)([ x(2,:) ; 4*9.8*sin(x(1, :)) - u ] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) - u]);

% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 3*x(1,:).^2.*x(2,:) - 1*x(2,:) - u] );
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
% liftFun = @(xx)( [xx; Encoder_VDP(xx)] - [zeros(2, 1); Encoder_VDP(zeros(2, 1))] );
% liftFun = @(xx)( [xx; xx(1).*xx(2); xx(1)^2.*xx(2); xx(2)^2.*xx(1)]);
% liftFun = @(xx)( [xx; xx(1).*xx(2); xx(1)^2.*xx(2); xx(2)^2.*xx(1); xx(2)^2.*xx(1).^2]);
liftFun = @(xx)( [xx]);
% liftFun = @(xx)( [xx;xx(1).*xx(2); xx(1)^2.*xx(2)]);
% liftFun = @(xx)( [xx; Hermite(1, xx); Hermite(2, xx)]);

Nsim = 100;
Ntraj = 100;

Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 2*rand([Nsim Ntraj]) - 1;
Xcurrent = (rand(n,Ntraj)*4 - 2);

N = Nsim * Ntraj / 2;

[A, B, C, Ylift, Xlift, X, U] = VDP_System_Generate(liftFun, 1, 2, 2, Ubig, Xcurrent, Nsim, Ntraj);

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


% Ylift3 = Ylift(:, 2*N + (1 : N));
% Xlift3 = Xlift(:, 2*N + (1 : N));
% X3 = X(:, 2*N + (1 : N));
% U3 = U(:, 2*N + (1 : N));
% W3 = [Ylift3 ; X3];
% V3 = [Xlift3; U3];
% V3Vt = V3*V3';
% W3Vt = W3*V3';
% M3 = W3Vt * pinv(V3Vt);

% A3 = M3(1:Nlift,1:Nlift);
% B3 = M3(1:Nlift,Nlift+1:end);
% C3 = M3(Nlift+1:end,1:Nlift);
% 
% 
% Ylift4 = Ylift(:, 3*N + (1 : N));
% Xlift4 = Xlift(:, 3*N + (1 : N));
% X4 = X(:, 3*N + (1 : N));
% U4 = U(:, 3*N + (1 : N));
% W4 = [Ylift4 ; X4];
% V4 = [Xlift4; U4];
% V4Vt = V4*V4';
% W4Vt = W4*V4';
% M4 = W4Vt * pinv(V4Vt);
% 
% A4 = M4(1:Nlift,1:Nlift);
% B4 = M4(1:Nlift,Nlift+1:end);
% C4 = M4(Nlift+1:end,1:Nlift);


% [A1, B1, C1, Ylift, Xlift, X, U] = VDP_System_Generate(liftFun, 1, 10, 2, Ubig, Xcurrent, Nsim, Ntraj);
% [A2, B2, C2, Ylift, Xlift, X, U] = VDP_System_Generate(liftFun, 1, 1, 2, Ubig, Xcurrent, Nsim, Ntraj);
% [A3, B3, C3, Ylift, Xlift, X, U] = VDP_System_Generate(liftFun, 0.1, 10, 2, Ubig, Xcurrent, Nsim, Ntraj);
% [A4, B4, C4, Ylift, Xlift, X, U] = VDP_System_Generate(liftFun, 0.1, 1, 2, Ubig, Xcurrent, Nsim, Ntraj);

% A2 = A1;
% B2 = B1;
% A3 = A1;
% B3 = B1;
% A4 = A1;
% B4 = B1;

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





%% Verification 
N = 1000;
x0_init = [-0.06; 0.83];
% x0_init = [-1.3; 1.5];
% x0_init = [0.2; 0.6];

% x0_init = 2 * rand(2, 1) - 1;
u = 2 * rand(1, N) - 1;
epsilon = 0.2;
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
    
    % A1 = 0.8 * rand(Nlift, Nlift);
    % A2 = 0.8 * rand(Nlift, Nlift);
    % N = 100;
    % NSpan = 1 : N;
    for i = 1 : N
        X1 = [X1 x1];
        X2 = [X2 x2];
        X = [X x0];
        X_Recon = [X_Recon x_rec];
    
        x_1 = C1 * A1 * x_rec + C1 *  B1 * u(i);
        x_2 = C2 * A2 * x_rec + C2 *  B2 * u(i);
%         x_3 = C3 * A3 * x_rec + C3 *  B3 * u(i);
%         x_4 = C4 * A4 * x_rec + C4 *  B4 * u(i);
        xlift = liftFun(x0);
    
        x0 = f_ud(0, x0, u(i));
        U0 = u(i);
        ylift = liftFun(x0);
    
    %     a_b = lsqlin([A1 * x_rec + B1 * u(i)  A2 * x_rec +  B2 * u(i)  A3 * x_rec +  B3 * u(i)  A4 * x_rec +  B4 * u(i) eye(Nlift, Nlift)], liftFun(x0),[ [-eye(4, 4) zeros(4, Nlift) ]; [eye(4, 4) zeros(4, Nlift)]; [ones(1, 4) zeros(1, Nlift)]], [zeros(4, 1); ones(4, 1); 1])
    %     a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1)
    %       a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ;  ones(1, 2)];], [zeros(2, 1);1])
        if j == 1
%             A1 = (1 + epsilon ) * A1;
%             B1 = (1 + epsilon ) * B1;
%             a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ;  ones(1, 2)];], [zeros(2, 1); 1 + epsilon])
            a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1)
%             a_b = lsqlin([x_1,  x_2,  x_3,  x_4], x0,[ [-eye(4, 4) ; ];], [zeros(4, 1);],  ones(1, 4), 1)

%             a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ;  ones(1, 2)];], [zeros(2, 1); 1])
        end
        if j == 2
%             a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1)
            a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ;  ones(1, 2)];], [zeros(2, 1); 1 ])
%             a_b = lsqlin([x_1,  x_2,  x_3,  x_4], x0,[ [-eye(4, 4) ; ones(1, 4)];], [zeros(4, 1);1])

%             a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1)

        end
    %     a_b = lsqlin([x_1,  x_2], x0,[ [-eye(2, 2) ; ];], [zeros(2, 1);],  ones(1, 2), 1 + epsilon)
    
    
    %     a_b = lsqlin([A1 * x_rec  A2 * x_rec  A3 * x_rec  A4 * x_rec   B1 * u(i)    B2 * u(i)    B3 * u(i)    B4 * u(i)], liftFun(x0), ...
    %                          [ [-eye(8, 8)  ]; [eye(8, 8) ]; ], [zeros(8, 1); ones(8, 1);], [[ones(1, 4) zeros(1, 4)]; [zeros(1, 4) ones(1, 4)]], [ 1; 1])
    
        a_b_Set = [a_b_Set a_b]; 
        
    %     rank([A1 * x_rec + B1 * u(i)  A2 * x_rec + B2 * u(i)  A3 * x_rec + B3 * u(i)  A4 * x_rec + B4 * u(i)])
    %     a_b = x0 * [A1 * x_rec + B1 * u(i); A2 * x_rec + B2 * u(i); A3 * x_rec + B3 * u(i); A4 * x_rec + B4 * u(i)]' * pinv([A1 * x_rec; A2 * x_rec; A3 * x_rec; A4 * x_rec; B1 * u(i);B2 * u(i);B3 * u(i);B4 * u(i)] * [A1 * x_rec; A2 * x_rec; A3 * x_rec; A4 * x_rec; B1 * u(i);B2 * u(i);B3 * u(i);B4 * u(i)]');
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
        
    %     x_rec = [a_b(1) * A1 * x_rec + a_b(5) * B1 * u(i) + a_b(2) * A2 * x_rec +  a_b(6) * B2 * u(i) +  a_b(3) * A3 * x_rec +  a_b(7) * B3 * u(i) +  a_b(4) * A4 * x_rec +  a_b(8) * B4 * u(i) ] ;
    %     x_rec = [A1 * x_rec  A2 * x_rec  A3 * x_rec  A4 * x_rec   B1 * u(i)    B2 * u(i)    B3 * u(i)    B4 * u(i)] * a_b;
    %     A = a_b(1) * A1 + a_b(2) * A2;
    %     B = a_b(1) * B1 + a_b(2) * B2;
%         x_rec = A1 * xlift + B1 * u(i);
        x_rec =  A * x_rec + B * u(i);
    
    
    end
    
    
    A1 = (1 + epsilon / 2) * A1;
    B1 = (1 + epsilon / 2) * B1;
    A2 = (1 + epsilon / 2) * A2;
    B2 = (1 + epsilon / 2) * B2;
%     A3 = diag(0.01 * ones(1, n));
%     B3 = zeros(n, 1);

%     [Ua, Sa, Va] = svd(A2);
%     [Ub, Sb, Vb] = svd(B2);
%     Sa = (1 + epsilon) * Sa;
%     Sb = (1 + epsilon) * Sb;
%     A2 = Ua*Sa*Va';
%     B2 = Ub*Sb*Vb';
    
%     [Ua, Sa, Va] = svd(A1);
%     [Ub, Sb, Vb] = svd(B1);
%     Sa = (1) * Sa;
%     Sb = (1) * Sb;
%     A1 = Ua*Sa*Va';
%     B1 = Ub*Sb*Vb';



    figure
    subplot 211
%     plot(NSpan, X1(1, :))
%     hold on
%     plot(NSpan, X2(1, :))
%     hold on
    plot(NSpan,X(1, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on
    plot(NSpan, X_Recon(1, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
    legend('Original System','LPV System')
    
    subplot 212 
%     plot(NSpan, X1(2, :))
%     hold on
%     plot(NSpan, X2(2, :))
%     hold on
    plot(NSpan,X(2, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
    hold on
    plot(NSpan, X_Recon(2, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
%     legend('Original System','LPV System')
    
    
    figure
    plot(NSpan, a_b_Set(1, :))
    hold on
    plot(NSpan, a_b_Set(2, :))
    legend('$p_1$', '$p_2$', 'Interpreter', 'latex')
    
    if j == 2
        figure
        plot(NSpan, Delta_p, 'LineWidth', 2, 'LineStyle', '-','Color','b')
        legend('$\Delta_p$', 'Interpreter', 'latex')
    end

    pause
    close all
end


%% Figure 5 Closed-loop responses for the time-varying system with input constraint; 

Control_Bound = 40;

Loops = 10;
X_Loop = [];
U_Loop = [];
for L = 1 : Loops

% A1 = A(0.1);
% A2 = A(10);
nx = size(A1, 1);
nu = size(B1, 2);
init = [-.5 + rand; -.5 + rand];
% init = [1; -1];
x0 = liftFun(init);
% Q1 = 1 * C1'*C1;
Q1 = diag([[10 10] zeros(1, Nlift - n)]);
R = 0.01;

F_Receding = -dlqr(A, B, Q1, R);
% 
% F_Receding = F_Receding_Static;
X_Receding = [init];
F_Set = [];
x_next = init;
U_Receding = [];
time = 200;
J_Set = [];
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

%     LMI4 = [Q                       (A3 * Q + B3 * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
%             A3 * Q + B3 * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
%             sqrt(Q1)*Q              zeros(nx, nx)        (gamma) * eye(nx, nx) zeros(nx, nu);
%             sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       (gamma) * eye(nu, nu)];
%     
%     LMI5 = [Q                       (A4 * Q + B4 * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
%             A4 * Q + B4 * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
%             sqrt(Q1)*Q              zeros(nx, nx)        (gamma) * eye(nx, nx) zeros(nx, nu);
%             sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       (gamma) * eye(nu, nu)];

% Input constraints  |uk| <= 2
% LMI0 = [X  Y;
%         Y' Q];
    LMI0 = [X  Y;
            Y' Q];

    ZEROS = 0;
    Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, Q >= ZEROS, LMI0 >= 0];
%     Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, Q >= ZEROS];


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
        
%     alpha = (10 - 0.1) * rand + 0.1;
%     A_Receding = A1;
    alpha = 0.2  + 0.8 * rand;
    beta = 8 + 2 * rand;
    gamma = 2 * rand;
%     f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) + u] );
%     f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 10 * x(1, :) - 2 * x(1, :).^3 + u] );
%     f_u =  @(t,x,u)([ x(2,:) ; -alpha*x(2, :) + beta * x(1, :) - 2 * x(1, :).^3 + u] );
%     f_u =  @(t,x,u)([ x(2,:) ;x(1, :).^2+ 0.15 .* u^3  + 0.1* (1 + x(2,:).^2*u + sin(0.1 .* u))] );

    k1 = @(t,x,u) ( f_u(t,x,u) );
    k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
    k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
    k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
    
    f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );
    %     x_next = A_Receding * x_next + B * F_Receding * x_next;
    Xlift1 = [Xlift1 liftFun(x_next)];
    Xlift2 = [Xlift2 liftFun(x_next)];
%     Xlift3 = [Xlift3 liftFun(x_next)];
%     Xlift4 = [Xlift4 liftFun(x_next)];
    J = J + x_next' * x_next + [F_Receding * liftFun(x_next)]'*[F_Receding * liftFun(x_next)];
    J_Set = [J_Set J];
    x_next = f_ud(0,x_next,F_Receding * liftFun(x_next))
%     x_next = f_ud(0,x_next,0)
    X_Receding = [X_Receding x_next];
    U_Receding = [U_Receding   F_Receding * liftFun(x_next)];
    F_Set  = [F_Set; F_Receding];
%     if norm(x_next) > 0.0001
%     
%         % Model update
%         X1 = [X1 x_next];
%         U1 = [U1 F_Receding * liftFun(x_next)];
%         Ylift1 = [Ylift1 liftFun(x_next)];
%         W1 = [Ylift1 ; X1];
%         V1 = [Xlift1; U1];
%         V1Vt = V1*V1';
%         W1Vt = W1*V1';
%         M1 = W1Vt * pinv(V1Vt);
%     
%         A1 = M1(1:Nlift,1:Nlift);
%         B1 = M1(1:Nlift,Nlift+1:end);
%         C1 = M1(Nlift+1:end,1:Nlift);
%     
%     
%         X2 = [X2 x_next];
%         U2 = [U2 F_Receding * liftFun(x_next)];
%         Ylift2 = [Ylift2 liftFun(x_next)];
%         W2 = [Ylift2 ; X2];
%         V2 = [Xlift2; U2];
%         V2Vt = V2*V2';
%         W2Vt = W2*V2';
%         M2 = W2Vt * pinv(V2Vt);
%     
%         A2 = M2(1:Nlift,1:Nlift);
%         B2 = M2(1:Nlift,Nlift+1:end);
%         C2 = M2(Nlift+1:end,1:Nlift);
% 
%         X3 = [X3 x_next];
%         U3 = [U3 F_Receding * liftFun(x_next)];
%         Ylift3 = [Ylift3 liftFun(x_next)];
%         W3 = [Ylift3 ; X3];
%         V3 = [Xlift3; U3];
%         V3Vt = V3*V3';
%         W3Vt = W3*V3';
%         M3 = W3Vt * pinv(V3Vt);
%     
%         A3 = M3(1:Nlift,1:Nlift);
%         B3 = M3(1:Nlift,Nlift+1:end);
%         C3 = M3(Nlift+1:end,1:Nlift);
% 
% 
%         X4 = [X4 x_next];
%         U4 = [U4 F_Receding * liftFun(x_next)];
%         Ylift4 = [Ylift4 liftFun(x_next)];
%         W4 = [Ylift4 ; X4];
%         V4 = [Xlift4; U4];
%         V4Vt = V4*V4';
%         W4Vt = W4*V4';
%         M4 = W4Vt * pinv(V4Vt);
%     
%         A4 = M4(1:Nlift,1:Nlift);
%         B4 = M4(1:Nlift,Nlift+1:end);
%         C4 = M4(Nlift+1:end,1:Nlift);
% 
%     end
end
    U_Loop = [U_Loop U_Receding];
    X_Loop = [X_Loop X_Receding];
end

% figure
% plot(0:0.1:time * 0.1,X_Receding(1,:), 'LineWidth', 3)
% hold on 
% plot(0:0.1:time * 0.1,X_Receding(2,:), 'LineWidth', 3)
% % axis([0 10,-0.2 1]);
% xlabel('$Time(sec)$','interpreter','latex');
% ylabel('$State$','interpreter','latex');
% title('$State\,\,trajectory$','interpreter','latex');



% figure
% plot(0:0.1:time * 0.1,U_Receding(1,:))
% % axis([0 10,-2 0.5]);
% xlabel('time(sec)');
% ylabel('$u$ (volts) ','interpreter','latex');
% title('Control signal $u$ (volts)','interpreter','latex');
Tspan = 0:0.1: time * 0.1;
NTime = size(Tspan, 2);
figure
for L = 0 : Loops - 1
    plot(Tspan,X_Loop(1, L * (time + 1) + 1 : (L + 1)* (time + 1)), 'LineWidth', 3, 'Color', 'b')
    hold on 
    plot(Tspan,X_Loop(2, L * (time + 1) + 1 : (L + 1)* (time + 1)), 'LineWidth', 3, 'Color', 'r')
end
% axis([0 10,-0.2 1]);
xlabel('$Time(sec)$','interpreter','latex');
ylabel('$State$','interpreter','latex');
title('$State\,\,trajectory$','interpreter','latex');


figure
for L = 0 : Loops - 1
    plot(Tspan(1: end - 1),U_Loop(1, L * (time + 1) + 1 : (L + 1)* (time)), 'LineWidth', 3)
    hold on 
end
plot(Tspan,Control_Bound * ones(1, NTime),'LineWidth', 3, 'LineStyle','--');
plot(Tspan,-Control_Bound* ones(1, NTime),'LineWidth', 3, 'LineStyle','--');

% axis([0 10,-2 0.5]);
xlabel('time(sec)');
ylabel('$u$ ','interpreter','latex');
title('Control signal $u$','interpreter','latex');

figure
plot(Tspan(1 : end - 1), J_Set,'LineWidth', 3, 'LineStyle','-')
xlabel('time(sec)');
ylabel('$J$ ','interpreter','latex');





