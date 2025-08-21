clc
clear
close all
% rng(66)
%% *************************** Dynamics ***********************************
% fd = lambda t, x, u: np.array([x[1, :], -0.5 * x[1, :] + x[0, :] - x[0, :] ** 3.0 + u])

% f_u =  @(t,x,u)([ 2*x(2,:) ; 2.0*x(2, :) - 2.0*x(1, :).^2.*x(2, :) - 0.8*x(1, :) + u] );
% f_u =  @(t,x,u)([ x(2,:) ; -0.7*x(2, :) + 2* x(1, :) - 2 *x(1, :).^3 + u] );
% f_u =  @(t,x,u)([ x(2,:) ; x(1, :) - x(2, :).^3 + u] );
% f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);


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



%% ************************** Collect data ********************************
% tic
disp('Starting data collection')
Nsim = 200;
Ntraj = 100;

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Random initial conditions
Xcurrent = (4*rand(n,Ntraj) - 2);

X = []; Y = []; U = [];
Xlift = []; Ylift = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end
% fprintf('Data collection DONE, time = %1.2f s \n', toc);

%% ************************** Basis functions *****************************

basisFunction = 'rbf';
% RBF centers
Nrbf = 2;
% cent = rand(n,Nrbf)*2 - 1;
[idx, cent_temp] = kmeans(X', Nrbf);
cent = cent_temp';
rbf_type = 'gauss'; 
% Lifting mapping - RBFs + the state itself
% liftFun = @(xx)( [xx;  rbf(xx,cent,rbf_type)] );
% Nlift = Nrbf + 2;
% liftFun = @(xx)( [xx; xx(1).*xx(2); xx(1)^2.*xx(2); xx(2)^2.*xx(1); xx(2)^2.*xx(1).^2]);
% liftFun = @(xx)( [xx; xx(1).*xx(2);  xx(1)^2.*xx(2);xx(2)^2.*xx(1);]);
% liftFun = @(xx)( [xx]);
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] - [zeros(2, 1);rbf(zeros(2, 1),cent,rbf_type)]);

% liftFun = @(x) Encoder_Duffing(x);
% Nlift = 8;
Nlift =size( liftFun(zeros(2, 1)), 1);




%% ******************************* Lift ***********************************
% 
% disp('Starting LIFTING')
% tic
% Xlift = liftFun(X);
% Ylift = liftFun(Y);
% fprintf('Lifting DONE, time = %1.2f s \n', toc);

for i = 1 : size(X, 2)
    Xlift = [Xlift liftFun(X(:, i))];
    Ylift = [Ylift liftFun(Y(:, i))];
end



%% ********************** Build predictor *********************************




Nsim = 200;
Ntraj = 100;

Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 2*rand([Nsim Ntraj]) - 1;
Xcurrent = (rand(n,Ntraj)*4 - 2);
[A1, B1, C1, Ylift1, Xlift1, X1, U1] = Duffing_System_Generate(liftFun, 1, 2, 2, Ubig, Xcurrent, Nsim, Ntraj);
% Ubig = 4*rand([Nsim Ntraj]) - 2;
% Xcurrent = (rand(n,Ntraj)*4 - 2);
% [A2, B2, C2, Ylift2, Xlift2, X2, U2] = Duffing_System_Generate(liftFun, 2, 2, 2,  Ubig, Xcurrent, Nsim, Ntraj);
% Ubig = 4*rand([Nsim Ntraj]) - 2;
% Xcurrent = (rand(n,Ntraj)*4 - 2);
% [A3, B3, C3, Ylift3, Xlift3, X3, U3] = Duffing_System_Generate(liftFun, 5, 2, 2,  Ubig, Xcurrent, Nsim, Ntraj);
% Ubig = 4*rand([Nsim Ntraj]) - 2;
% Xcurrent = (rand(n,Ntraj)*4 - 2);
% [A4, B4, C4, Ylift4, Xlift4, X4, U4] = Duffing_System_Generate(liftFun, 7, 2, 2,  Ubig, Xcurrent, Nsim, Ntraj);
% 
% 
% N = Nsim;
% % NSpan = 1 * deltaT : deltaT : Nsim * deltaT;
% NSpan = 1 : Nsim;
% 
% x0 = [0.5; 1.2];
% % u = rand(1, N);
% u = prbs(15, Nsim);
% X1 = [];
% X2 = [];
% X3 = [];
% X4 = [];
% X = [];
% x1 = liftFun(x0);
% x2 = liftFun(x0);
% x3 = liftFun(x0);
% x4 = liftFun(x0);
% for i = 1 : N
%     X1 = [X1 x1];
%     X2 = [X2 x2];
%     X3 = [X3 x3];
%     X4 = [X4 x4];
%     X = [X x0];
% 
%     x0 = f_ud(0, x0, u(i));
%     x1 = A1 * x1 + B1 * u(i);
%     x2 = A2 * x2 + B2 * u(i);
%     x3 = A3 * x3 + B3 * u(i);
%     x4 = A4 * x4 + B4 * u(i);
% end
% 
% figure
% plot(NSpan, X1(1, :))
% hold on
% plot(NSpan, X2(1, :))
% hold on
% plot(NSpan, X3(1, :))
% hold on
% plot(NSpan, X4(1, :))
% hold on
% plot(NSpan, X(1, :))
% 
% legend('A1', 'A2', 'A3', 'A4', 'Origin')
% 
% 
% figure
% plot(NSpan, X1(2, :))
% hold on
% plot(NSpan, X2(2, :))
% hold on
% plot(NSpan, X3(2, :))
% hold on
% plot(NSpan, X4(2, :))
% hold on
% plot(NSpan, X(2, :))
% 
% 
% legend('A1', 'A2', 'A3', 'A4', 'Origin')
% 
% 
% figure
% plot(X1(1, :), X1(2, :))
% hold on
% plot(X2(1, :),X2(2, :))
% hold on
% plot(X3(1, :),X3(2, :))
% hold on
% plot(X4(1, :),X4(2, :))
% hold on
% plot(X(1, :),X(2, :))
% 
% 
% legend('A1', 'A2', 'A3', 'A4', 'Origin')
% 
% Theta1 = [];
% Theta2 = [];
% for i = 1 : size(X1, 2)
% %     theta1 = lsqlin([[C1 * X1(:, i)]' [C2 *X2(:, i)]' [C3 *X3(:, i)]' [C4 *X4(:, i)]'], X(:, i)', [-eye(8); eye(8)], [zeros(4, 1); 0.5 * ones(4, 1)])
% %     theta1 = lsqlin([X1(1, i)' X2(1, i)' X3(1, i)' X4(1, i)'], X(1, i)', [-eye(4); eye(4); ones(1, 4)], [zeros(4, 1); 1 * ones(4, 1); 1])
%     theta1 = lsqlin([X1(1, i)' X2(1, i)' X3(1, i)' X4(1, i)'], X(1, i)', [-eye(4)], [zeros(4, 1)], ones(1, 4), 1)
%     theta2 = lsqlin([X1(2, i)' X2(2, i)' X3(2, i)' X4(2, i)'], X(2, i)', [-eye(4)], [zeros(4, 1)], ones(1, 4), 1)
%     Theta1 = [Theta1 theta1];
%     Theta2 = [Theta2 theta2];
% end
% % theta1 = lsqlin([X1(1:2, i)' X2(1:2, i)' X3(1:2, i)' X4(1:2, i)'], X(1:2, i)', [-eye(4); eye(4); ones(1, 4)], [zeros(4, 1); 1 * ones(4, 1); 1])
% 
% % theta = [X] * pinv([[C1 * X1] ;[C2 *X2] ;[C3 *X3] ;[C4 *X4]])
% p = sum(Theta1, 2)  / size(Theta1, 2);
% 
% X1_Recon = [];
% for i = 1 : size(X1, 2)
%     X1_Recon = [X1_Recon X1(1, i)' * Theta1(1, i) +  X2(1, i)' * Theta1(2, i) +  X3(1, i)'  * Theta1(3, i) + X4(1, i)' * Theta1(4, i)];
% end
% 
% X2_Recon = [];
% for i = 1 : size(X1, 2)
%     X2_Recon = [X2_Recon X1(2, i)' * Theta2(1, i) +  X2(2, i)' * Theta2(2, i) +  X3(2, i)'  * Theta2(3, i) + X4(2, i)' * Theta2(4, i)];
% end
% 
% figure
% plot(NSpan, X(1, :))
% legend('Origin')
% 
% hold on
% plot(NSpan, X1_Recon(1, :))
% legend('X1 Reconstruction')
% 
% 
% figure
% plot(NSpan, X(2, :))
% legend('Origin')
% 
% hold on
% plot(NSpan, X2_Recon(1, :))
% legend('X2 Reconstruction')
% 
% 
% N = 3000;
% NSpan = 1 : N;
% x0 = [0.5; 1.2];
% u = 2 * rand(1, N) - 1;
% % u = prbs(15, N);
% X1 = [];
% X2 = [];
% X3 = [];
% X4 = [];
% X = [];
% X_Recon = [];
% x1 = liftFun(x0);
% x2 = liftFun(x0);
% x3 = liftFun(x0);
% x4 = liftFun(x0);
% x_rec = liftFun(x0);
% a_b_Set = [];
% for i = 1 : N
%     X1 = [X1 x1];
%     X2 = [X2 x2];
%     X3 = [X3 x3];
%     X4 = [X4 x4];
%     X = [X x0];
%     X_Recon = [X_Recon x_rec];
% 
%     x0 = f_ud(0, x0, u(i));
% %     a_b = lsqlin([A1 * x_rec + B1 * u(i)  A2 * x_rec +  B2 * u(i)  A3 * x_rec +  B3 * u(i)  A4 * x_rec +  B4 * u(i) eye(Nlift, Nlift)], liftFun(x0),[ [-eye(4, 4) zeros(4, Nlift) ]; [eye(4, 4) zeros(4, Nlift)]; [ones(1, 4) zeros(1, Nlift)]], [zeros(4, 1); ones(4, 1); 1])
% %     a_b = lsqlin([A1 * x_rec + B1 * u(i)  A2 * x_rec +  B2 * u(i)  A3 * x_rec +  B3 * u(i)  A4 * x_rec +  B4 * u(i)], liftFun(x0),[ [-eye(4, 4)  ]; [eye(4, 4) ]], [zeros(4, 1); ones(4, 1)])
% %     a_b = lsqlin([A1 * x_rec  A2 * x_rec  A3 * x_rec  A4 * x_rec   B1 * u(i)    B2 * u(i)    B3 * u(i)    B4 * u(i)], liftFun(x0), ...
% %                          [ [-eye(8, 8)  ]; [eye(8, 8) ]; ], [zeros(8, 1); ones(8, 1);], [[ones(1, 4) zeros(1, 4)]; [zeros(1, 4) ones(1, 4)]], [ 1; 1])
%     a_b = lsqlin([A1 * x_rec + B1 * u(i)  A2 * x_rec +  B2 * u(i)  A3 * x_rec +  B3 * u(i)  A4 * x_rec +  B4 * u(i)], liftFun(x0),[ [-eye(4, 4)  ]; [eye(4, 4) ]], [zeros(4, 1); ones(4, 1)])
% 
%     a_b_Set = [a_b_Set a_b]; 
% %     rank([A1 * x_rec + B1 * u(i)  A2 * x_rec + B2 * u(i)  A3 * x_rec + B3 * u(i)  A4 * x_rec + B4 * u(i)])
% %     a_b = x0 * [A1 * x_rec + B1 * u(i); A2 * x_rec + B2 * u(i); A3 * x_rec + B3 * u(i); A4 * x_rec + B4 * u(i)]' * pinv([A1 * x_rec; A2 * x_rec; A3 * x_rec; A4 * x_rec; B1 * u(i);B2 * u(i);B3 * u(i);B4 * u(i)] * [A1 * x_rec; A2 * x_rec; A3 * x_rec; A4 * x_rec; B1 * u(i);B2 * u(i);B3 * u(i);B4 * u(i)]');
%     x1 = A1 * x1 + B1 * u(i);
%     x2 = A2 * x2 + B2 * u(i);
%     x3 = A3 * x3 + B3 * u(i);
%     x4 = A4 * x4 + B4 * u(i);
% %     x_rec = [a_b(1) * A1 * x_rec + a_b(5) * B1 * u(i) + a_b(2) * A2 * x_rec +  a_b(6) * B2 * u(i) +  a_b(3) * A3 * x_rec +  a_b(7) * B3 * u(i) +  a_b(4) * A4 * x_rec +  a_b(8) * B4 * u(i) ] ;
% %     x_rec = [A1 * x_rec  A2 * x_rec  A3 * x_rec  A4 * x_rec   B1 * u(i)    B2 * u(i)    B3 * u(i)    B4 * u(i)] * a_b;
%     x_rec = [a_b(1) * A1 * x_rec + a_b(1) * B1 * u(i) + a_b(2) * A2 * x_rec +  a_b(2) * B2 * u(i) +  a_b(3) * A3 * x_rec +  a_b(3) * B3 * u(i) +  a_b(4) * A4 * x_rec +  a_b(4) * B4 * u(i) ] ;
% end
% 
% figure
% plot(NSpan, X1(1, :))
% hold on
% plot(NSpan, X2(1, :))
% hold on
% plot(NSpan, X3(1, :))
% hold on
% plot(NSpan, X4(1, :))
% hold on
% plot(NSpan,X(1, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
% hold on
% plot(NSpan, X_Recon(1, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
% % A = p(1) * A1 + p(2) * A2+ p(3) * A3 + p(4) * A4;
% % B = p(1) * B1 + p(2) * B2+ p(3) * B3 + p(4) * B4;
% legend('A1', 'A2', 'A3', 'A4', 'Origin','Recon')
% 
% 
% figure
% plot(NSpan, X1(2, :))
% hold on
% plot(NSpan, X2(2, :))
% hold on
% plot(NSpan, X3(2, :))
% hold on
% plot(NSpan, X4(2, :))
% hold on
% plot(NSpan, X(2, :), 'LineWidth', 3, 'LineStyle', '--','Color','r')
% hold on
% plot(NSpan, X_Recon(2, :), 'LineWidth', 2, 'LineStyle', '-','Color','b')
% % A = p(1) * A1 + p(2) * A2+ p(3) * A3 + p(4) * A4;
% % B = p(1) * B1 + p(2) * B2+ p(3) * B3 + p(4) * B4;
% legend('A1', 'A2', 'A3', 'A4', 'Origin','Recon')


A = A1;
B = B1;
C = C1;

% disp('Starting REGRESSION')
% tic
% W = [Ylift ; X];
% V = [Xlift; U];
% VVt = V*V';
% WVt = W*V';
% M = WVt * pinv(VVt); % Matrix [A B; C 0]
% A = M(1:Nlift,1:Nlift);
% B = M(1:Nlift,Nlift+1:end);
% C = M(Nlift+1:end,1:Nlift);
% 
% fprintf('Regression done, time = %1.2f s \n', toc);


nx = size(B, 1);
nu = size(B, 2);
Cy = eye(n);
%     Cy = [eye(size(Yr, 1))];
% N = 10;
N = 30;

Q = diag([10 10]);
R = 0.001;
Yr = [0; 0];
Q_bar = kron(eye(N), Q);
R_bar = kron(eye(N), R);
% Q_bar = kron(eye(N), Qlift);
% R_bar = kron(eye(N), Rlift);
% Q_bar(end - Nlift + 1 : end, end - Nlift + 1 : end) = P;
x0 = [1; -1];
Lift_xu = [liftFun(x0)];
Shift_Matrix = kron([zeros(1, N); [eye(N - 1) zeros(N - 1, 1)]], eye(nu));
Compact_Form1 = [];
for i = 1 : N
    Compact_Form1 = [Compact_Form1; Cy * C * A^i];
end
Compact_Form2 = [];
for i = 1 : N
    vector_Temp = [];
    for j = 1 : N
        vector_Temp = [Cy * C * A^(j - 1)*B  vector_Temp];
    end
    Compact_Form2 = [vector_Temp * Shift_Matrix^(i - 1); Compact_Form2];
end
Yr = kron(ones(N, 1), Yr);
p = size(Yr, 1);
H = (Compact_Form2)' * Q_bar * (Compact_Form2) + R_bar;
H = (H+H')/2;
f = 2 .* (Compact_Form1 * Lift_xu)' * Q_bar * Compact_Form2 - 2 .* Yr' * Q_bar * Compact_Form2;
[U0_Set, feval] = quadprog(2.*H, f, [], [], [], [], -2*ones(p, 1), 2*ones(p, 1));
U0 = U0_Set(1 : nu)

X_Collection = [];
U_Collection = [];
Steps = 200;
J = 0;
J_Set = [];
Ref_Plot = [];
time = 0;
x_max = kron(ones(N , 1), [1.2; 1.0]);
x_min = kron(ones(N , 1), [-1.0; -1.0]);
AIneq = Compact_Form2;
AIneq = [AIneq; -Compact_Form2];
bIneq = [[x_max - Compact_Form1 * Lift_xu];
              [-x_min + Compact_Form1 * Lift_xu]];
for i = 1 : Steps
    %% MPC Solve
    AIneq = Compact_Form2;
    AIneq = [AIneq; -Compact_Form2];
    bIneq = [[x_max - Compact_Form1 * Lift_xu];
    [-x_min + Compact_Form1 * Lift_xu]];
    if mod(i, 1) == 0
        flag = 0;
        f = 2 .* (Compact_Form1 * Lift_xu)' * Q_bar * Compact_Form2 - 2 .* Yr' * Q_bar * Compact_Form2;
        [U0_Set, fval] = quadprog(2.*H, f, AIneq, bIneq, [], [], -4*ones(N, 1), 4*ones(N, 1));
        feval = fval + (Compact_Form1 * Lift_xu)' * Q_bar * (Compact_Form1 * Lift_xu)
    end
    Ref_Plot = [Ref_Plot; 0];
    U0 = U0_Set(flag * nu + 1 : flag * nu + nu);

    if i > 10
%         f_u =  @(t,x,u)([ 2*x(2,:) ; 2.0*x(2, :) - 2.0*x(1, :).^2.*x(2, :) - 0.8*x(1, :) + u] );
%          f_u =  @(t,x,u)([ 1.0*x(2,:) ; 5.0*x(2, :) - 2.0*x(1, :).^2.*x(2, :) - 1.2*x(1, :) + u] );
        % f = lambda t, x, u: np.array([x[1, :], -10.0 * 0.5 * x[1, :] + 2.0 * x[0, :] - 0.5 * x[0, :] ** 3.0 + u])
%         f_u =  @(t,x,u)([ x(2,:) ; - 10*x(2,:)*0.5 + 2.0 * x(1, :) - 0.5 * x(1, :).^3  + u] );
        k1 = @(t,x,u) (  f_u(t,x,u) );
        k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
        k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
        k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
        f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );
    end

    X = [X x0];
    xlift = [liftFun([x0])];
    Xlift = [Xlift [liftFun([x0])]];
    X_Collection = [X_Collection x0];
    x = x0;
    x0 = f_ud(0, x0, U0)
    y = x0;

    ylift = [liftFun([x0])];
    Ylift = [Ylift [liftFun([x0])]];

    U = [U U0];
    U_Collection = [U_Collection U0];
    tic
    lambda = 1.0;

%     a_b = lsqlin([A1 * xlift + B1 * U0   A2 * xlift +  B2 * U0  A3 * xlift +  B3 * U0  A4 * xlift +  B4 * U0], liftFun(x0),[ [-eye(4, 4)  ]; [eye(4, 4) ]], [zeros(4, 1); ones(4, 1)])
%     
%     A =  a_b(1) * A1 + a_b(2) * A2 +  a_b(3) * A3  +  a_b(4) * A4 ;
%     B =  a_b(1) * B1 + a_b(2) * B2 +  a_b(3) * B3  +  a_b(4) * B4 ;


%  %自定义满秩
%     if i == 1 
%         K_A = zeros(Nlift, Nlift + nu);
% %         K_A = V;
%         invK_G = 1e4 * eye(Nlift + nu);
%  %        invK_G = pinv(invK_G);
%  %        invK_G = pinv(W);
%         invK_G = 1 / lambda * invK_G - 1 / lambda * (invK_G * [xlift; U0] * [xlift; U0]'* invK_G) / (lambda + [xlift; U0]' * invK_G * [xlift; U0]);
%        K_A = K_A + [ylift] * [xlift; U0]';
%     else
%         invK_G = 1 / lambda * invK_G - 1 / lambda * (invK_G * [xlift; U0] * [xlift; U0]'* invK_G) / (lambda + [xlift; U0]' * invK_G * [xlift; U0]);
%         K_A = K_A + [ylift] * [xlift; U0]';
%     end
%     Kext = K_A * invK_G;
%     A = Kext(:, 1:Nlift);
%     B = Kext(:, Nlift + 1:end);

    
%     if i == 1 
%         bar_X = zeros(n, Nlift);
% %         K_A = V;
%         bar_X = bar_X + x * xlift';
%         bar_Q = 10000 * eye(Nlift);
%         bar_Q = 1 / lambda * bar_Q - 1 / lambda * (bar_Q * xlift * xlift' * bar_Q) / (lambda + xlift' * bar_Q * xlift);
% 
%     else
%         bar_X = bar_X + x * liftFun([x])';
%         bar_Q = 1 / lambda * bar_Q - 1 / lambda * (bar_Q * xlift * xlift' * bar_Q) / (lambda + xlift' * bar_Q * xlift);
%     end
%     C = bar_X * bar_Q;
%     t2 = toc
%     time = t2 + time;
%     % 储存法
%     tic
% 
%     W = [Ylift];
%     V = [Xlift; U];
%     VVt = V*V';
%     WVt = W*V';
%     M = WVt * pinv(VVt); % Matrix [A B; C 0]
%     A = M(1:Nlift,1:Nlift);
%     B = M(1:Nlift,Nlift+1:end);
%     C = (X * Xlift' ) * pinv(Xlift * Xlift');
%     t2 = toc
%     time = time + t2;
    J_Set = [J_Set J];
    J = J + x0'*x0 + U0'*U0;
    Shift_Matrix = kron([zeros(1, N); [eye(N - 1) zeros(N - 1, 1)]], eye(nu));
    Compact_Form1 = [];
    for k = 1 : N
        Compact_Form1 = [Compact_Form1; Cy * C * A^k];
    end
    Compact_Form2 = [];
    for k = 1 : N
        vector_Temp = [];
        for j = 1 : N
            vector_Temp = [Cy * C * A^(j - 1)*B  vector_Temp];
        end
        Compact_Form2 = [vector_Temp * Shift_Matrix^(k - 1); Compact_Form2];
    end
    p = size(Yr, 1);
    H = (Compact_Form2)' * Q_bar * (Compact_Form2) + R_bar;
    H = (H+H')/2;
    f = 2 .* (Compact_Form1 * Lift_xu)' * Q_bar * Compact_Form2 - 2 .* Yr' * Q_bar * Compact_Form2;
    
    Lift_xu = liftFun(x0);
end

% MSE = (X_Collection(1, :) - Ref_Plot') * (X_Collection(1, :) - Ref_Plot')' / Steps
Steady_Error = abs(X_Collection(1, end) - Ref_Plot(1))

figure 
plot(X_Collection(1, :),'LineStyle','-','LineWidth',2.0)
hold on
plot(X_Collection(2, :),'LineStyle','-','LineWidth',2.0)
% plot(x_max(1) * ones(size(X_Collection, 2), 1),'LineStyle','--','LineWidth',2.0)

figure 
plot(U_Collection(1, :),'LineStyle','-','LineWidth',2.0)


figure 
plot(J_Set(1, :),'LineStyle','-','LineWidth',2.0)







