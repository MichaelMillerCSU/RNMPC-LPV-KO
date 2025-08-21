clc
clear
close all

r = randi(5000);
n = 2;
Nd = 2;
% NK = n;
% rng(2958) % Nd = 10 very good
% rng(318) % All effective
Taylor = 0;
rng(1339) 


Det = [];
Eig = [];
Eig_max = [];
%% System set-up
%% Rotating arm 
J = 0.0292;
cm = 16;
Tc = 0.416;
Ts = 0.4657;
vs = 0.2;
sigma2 = 0.0135;
% f_u =  @(t,x,u)([ x(2,:) ; -Tc / J .* sign(x(2,:)) - (Ts - Tc) ./ J .* exp(-(x(2,:) / vs)) .* sign(x(2,:)) - sigma2 ./ J .* x(2,:) + cm / J .* u]);

%% Inverted pendulum on a cart
g = 9.8;
m = 2;
M = 8;
a = 1 / (m + M);
l = 0.5;

% f_u =  @(t,x,u)([ x(2,:) ; (g * sin(x(1, :)) - 0.5 * a * m * l * x(2,:).^2 .* sin(2 * x(1,:)) ) / (4 * l / 3 - a * m * l * cos(x(1, :).^2)) - a * cos(x(1,:)) ./ (4 * l / 3 - a * m * l * cos(x(1, :).^2)) .* u]);
f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2* x(1, :) - 2 *x(1, :).^3 + u] );

% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)([ 1*x(1,:) ; 4 * 9.8 * sin(x(1,:)) - 3 * u.* cos(x(1, :))] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 10*x(1,:).^2.*x(2,:) - 3*x(2,:) + u] );
deltaT = 0.05;
liftFun = @(xx)( [xx]);
% liftFun = @(xx)( [xx; xx(1).*xx(2)]);

basisFunction = 'rbf';
Nrbf = 2;
% cent = rand(2,Nrbf)*4 - 2;
cent = [[0.6486; -1.7319], [-0.5527;1.7656]];
rbf_type = 'gauss';
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)]- [zeros(2, 1); rbf(zeros(2, 1), cent,rbf_type)]);


%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

Nsim = 200;
Ntraj = 100;
Nlift =size( liftFun(zeros(2, 1)), 1);
Ubig = 2*rand([Nsim Ntraj]) - 1;
Xcurrent = (rand(n,Ntraj)*4 - 2);
[Alift, Blift, Clift, Ylift, Xlift, X, Y, U] = System_Data_Obtain(liftFun, f_u, deltaT, Ubig, Xcurrent, Nsim, Ntraj);
N = Nsim * Ntraj / 2; 

% A_Save = A;
% B_Save = B;
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
% A{1} = Alift;
% B{1} = Blift;

A{1} = A1;
B{1} = B1;

A{2} = 1.01 * A2;
B{2} = 1.01 * B2;

A_Rec = Alift;
B_Rec = Blift;
C = [1 0];
A{end + 1} = [0 * rand(Nlift, Nlift)];
B{end + 1} = 0 * rand(Nlift, 1);
Nd = Nd + 1; 

N = 1000;
Tspan = (1 : N) * deltaT;

%% Testing the controller on initial point
% load Coefficient_1.mat
% P_i_Set = P_i;
K_lqr = [-17.7137  -20.4804]; % LQR control for nonlinear system with VDP (0.8)
Times = 100;
check = 1;
E = [];
EigEA = [];
EigEB = [];
Err = [];
u = 1;
RankAy = [];
RankABy = [];
by_Set_No = [];
by_Set = [];
lambda = 0;
Convex_hull_data = [];
X_all = [];
U_all = [];
X_test = [];
Scheduling_Set = [];

for k = 1 : Times
    x = liftFun(2 * rand(2, 1) - 1);
%     x = liftFun([-1; 1]);

    x_rec = x;
    X_test = [x];
    X_recon_test = [x_rec];
    P_i = [];
    U = [];
    u = 0;
    J = 0 ;
    for i = 1 : N
        w = 0;
        u = 2 * rand - 1;
%         [~, u, ~] = Cell_LMI_Solve(A_Rec, B_Rec, A, B, x, Nd)
%         [~, u, ~] = Cell_LMI_Solve(A{2}, B{2}, A, B, x, Nd)

%         w = 0.001 * rand - 0.001 / 2;
%         u = K_lqr * x;
%         u = sin(0.1 * i);
%         if k > 20
%             u = K_lqr * x;
%         else
%             u = 4 * rand - 2;
%         end
        for j = 1 : Nd
            x_candidate{j} = A{j} * x_rec + B{j} * u;
            A_vec_candidation{j} = [vec([A{j}]); vec([B{j}])];
        end
        x_next = f_ud(0, Clift * x, u) + w;
        J = J + x_next'*x_next + u'*u
        x = liftFun(x_next)
        
        M = cell2mat(x_candidate);
        MA = [cell2mat(A_vec_candidation)];
        p_i = lsqlin(M, x,[ [-eye(Nd, Nd) ; ];], [zeros(Nd, 1);],  ones(1, Nd), 1);
%         p_i = lsqlin(M, x,[ [-eye(Nd, Nd) ; ones(1, Nd)];], [zeros(Nd, 1); 1]);

        U = [U u];

        sum(p_i)

        P_i = [P_i p_i];
        A_Nd = cell2mat(A);
        B_Nd = cell2mat(B);
        % Taylor expansion
        if Taylor == 1
            R_Nd = cell2mat(residual);
            R_Rec = R_Nd * p_i;
        end
    
        A_temp = A_Nd * kron(p_i, eye(Nlift));
        B_temp = B_Nd * p_i;
        
% 
        A_Rec = A_temp;
        B_Rec = B_temp;

        A_Rec_Set{i} = A_Rec;
        B_Rec_Set{i} = B_Rec;

%         [~, u, ~] = Cell_LMI_Solve(A_Rec, B_Rec, A, B, x, Nd)
        x_rec = A_Rec * x_rec + B_Rec * u;
        

%         Now_e = liftFun(x) - x_rec

        % Taylor expansion
        if Taylor == 1
            x_rec = A_Rec * x_rec + B_Rec * u + R_Rec;
        end
        X_test = [X_test x];
        X_recon_test = [X_recon_test x_rec];

        disp(i / N * 100)
    end

    Error = norm(X_test - X_recon_test) / norm(X_test)
    Err = [Err Error];

%     figure
%     subplot(311)
%     plot(X_test(1, :))
% 
%     hold on 
%     plot(X_recon_test(1, :))
%     legend('Original', 'LPV')
%     ylabel('x_1')
% 
% 
% 
%     subplot(312)
%     plot(X_test(2, :))
%     legend('Original')
% 
%     hold on 
%  
%     plot(X_recon_test(2, :))
%     legend('Original', 'LPV')
%     ylabel('x_2')
% 
%     subplot(313)
% 
%     mesh(Tspan, 1 : Nd, P_i)
%     xlabel('Time(sec)')
%     ylabel('i-th vertex')
%     zlabel('p_i')
% %     ylabel('x_2')
% 
% % %     figure
% % %     plot(U(:))
% % %     legend('Control input')
%     pause(1)
    X_all = [X_all X_test(:, 1 : end - 1)];
    U_all = [U_all U];
    Scheduling_Set = [Scheduling_Set P_i];
end

% pause()
%% Part of NN Learning Scheduling
close all
x_num = size(X_all, 2);

%% 设置训练数据和预测数据
% input_train = [X_all; U_all];
input_train = [X_all];

%% 训练样本数据归一化
[inputn,inputps]=mapminmax(input_train);%归一化到[-1,1]之间，inputps用来作下一次同样的归一化

% inputn = input_train;
% outputn = output_train;

% net1=newff(inputn,outputn,[20 20 20],{'tansig','purelin'},'traingdx');% 建立模型，传递函数使用purelin，采用梯度下降法训练
Net_Vertices = {};
% parpool("local",12);

for i = 1 : 1
    output_train = Scheduling_Set;
    % 构建BP神经网络
    Net_Vertices{i} = feedforwardnet([20 Nd], 'trainbr');
    Net_Vertices{i}.layers{1}.transferFcn = 'poslin'; % 
    Net_Vertices{i}.layers{2}.transferFcn = 'softmax'; % 

    % 网络参数配置（ 训练次数，学习速率，训练目标最小误差等）
    Net_Vertices{i}.trainParam.epochs=1e3;         % 训练次数，这里设置为1000次
    Net_Vertices{i}.trainParam.lr=0.001;                   % 学习速率，这里设置为0.01
    Net_Vertices{i}.trainParam.goal=1e-3;                    % 训练目标最小误差，这里设置为0.00001
%     lambda = 0.001;
    Net_Vertices{i}.performFcn = 'mse'; 
    Net_Vertices{i}.performParam.regularization = 1e-7; 
    Net_Vertices{i}.performParam.normalization = 'standard';
    Net_Vertices{i}=train(Net_Vertices{i},input_train,output_train);%开始训练，其中inputn,outputn分别为输入输出样本

%     Net_Vertices{i} = newrb(input_train, output_train, 1e-5);
end


% Tar=sim(net1,inputn_test); %用训练好的模型进行仿真
%% Simulation using NN
load Duffing_New.mat
close all
N = 200;
Tspan = (1 : N) * deltaT;
A_Rec = A{1};
B_Rec = B{1};
for k = 1 : 1
%     x = [-1 + 0.2 * k; 1 - 0.2 * k];
%     x = 2 * rand(2, 1) - 1;
    x = liftFun([-1; 0.5]);
%     x = [-1.5; 1];
    x_rec = x;
    X_test = [x];
    X_recon_test = [x_rec];
    P_i = [];
    U = [];
    J = 0;
    J_Set = [];
    for i = 1 : N
        w = 0;
%         w = 0.001 * rand - 0.001 / 2;
%         u = [-1 -2] * x;
%         u = 2 * rand - 1;
%         u = sin(0.1 * i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         p_i = sim(Net_Vertices{1},[x]);
%         A_temp = A_Nd * kron(p_i, eye(Nlift));
%         B_temp = B_Nd * p_i;
%         
%         A_Rec = A_temp;
%         B_Rec = B_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, u, ~] = Cell_LMI_Solve(A_Rec, B_Rec, A, B, x, Nd)
%         [~, u_rec, ~] = Cell_LMI_Solve(A_Rec, B_Rec, A, B, x_rec, Nd)

        for j = 1 : Nd
            x_candidate{j} = A{j} * x + B{j} * u;
            A_vec_candidation{j} = [vec([A{j}]); vec([B{j}])];
        end

%         p_i = [];
%         for j = 1 : Nd
%             p_i = [p_i; sim(Net_Vertices{j},[x; u])];
%         end
%         p_i = p_i ./ min(p_i) ;
%         p_i = p_i ./sum(p_i);

        x_next = f_ud(0, Clift * x, u) + w;
        J = J + x_next'* x_next + u' *u
        J_Set = [J_Set J]
        x = liftFun(x_next)

        M = cell2mat(x_candidate);
        MA = [cell2mat(A_vec_candidation)];
        p_i = lsqlin(M, x,[ [-eye(Nd, Nd) ; ];], [zeros(Nd, 1);],  ones(1, Nd), 1);
        
        U = [U u];

        sum_p = sum(p_i)

        
        P_i = [P_i p_i];
        A_Nd = cell2mat(A);
        B_Nd = cell2mat(B);
        if Taylor == 1
            R_Nd = cell2mat(residual);
            R_Rec = R_Nd * p_i;
        end
    
        A_temp = A_Nd * kron(p_i, eye(Nlift));
        B_temp = B_Nd * p_i;
        
        A_Rec = A_temp;
        B_Rec = B_temp;


        A_Rec_Set{i} = A_Rec;
        B_Rec_Set{i} = B_Rec;
        x_rec = A_Rec * x + B_Rec * u;
%         Now_e = x - x_rec

        % Taylor expansion
        if Taylor == 1
            x_rec = A_Rec * x_rec + B_Rec * u + R_Rec;
        end
        X_test = [X_test x];
        X_recon_test = [X_recon_test x_rec];

        disp(i)
    end

    Error = norm(X_test - X_recon_test, 'fro') / norm(X_test, 'fro')
    Err = [Err Error];


    
    figure
    subplot(311)
    plot(X_test(1, :))

    hold on 
    plot(X_recon_test(1, :))
    legend('Original', 'LPV')
    ylabel('x_1')

    subplot(312)
    plot(X_test(2, :))
    legend('Original')

    hold on 
    plot(X_recon_test(2, :))
    legend('Original', 'LPV')
    ylabel('x_2')

    subplot(313)

    mesh(Tspan, 1 : Nd, P_i)
    xlabel('Time(sec)')
    ylabel('i-th vertex')
    zlabel('p_i')

    figure
    plot(U(:))
    legend('Control input')

end


% %%  NN Quasi min-max MPC
% close all
% N = 200;
% J = 0;
% Tspan = (1 : N) * deltaT;
% A_Rec = A{1};
% B_Rec = B{1};
% for k = 1 : 1
% %     x = [-1 + 0.2 * k; 1 - 0.2 * k];
% %     x = 4 * rand(2, 1) - 2;
%     x = [1; -0.5];
%     x_rec = x;
%     X_test = [x];
%     X_recon_test = [x_rec];
%     P_i = [];
%     U = [];
%     for i = 1 : N
%         w = 0;
% %         w = 0.001 * rand - 0.001 / 2;
% %         u = [-1 -2] * x;
% %         u = 2 * rand - 1;
% %         u = sin(0.1 * i);
%         [~, u, ~] = Cell_LMI_Solve(A_Rec, B_Rec, A, B, x, Nd)
% 
%         p_i = sim(Net_Vertices{1},[x; u]);
% 
%         x_next = f_ud(0, Clift * x, u) + w;
%         x = liftFun(x_next)
% 
%         M = cell2mat(x_candidate);
%         MA = [cell2mat(A_vec_candidation)];
% %         p_i = lsqlin(M, x,[ [-eye(Nd, Nd) ; ];], [zeros(Nd, 1);],  ones(1, Nd), 1);
%         
%         U = [U u];
% 
%         sum_p = sum(p_i)
% 
%         
%         P_i = [P_i p_i];
%         A_Nd = cell2mat(A);
%         B_Nd = cell2mat(B);
%         if Taylor == 1
%             R_Nd = cell2mat(residual);
%             R_Rec = R_Nd * p_i;
%         end
%         J = J + x'* x + u' *u
% 
%         A_temp = A_Nd * kron(p_i, eye(Nlift));
%         B_temp = B_Nd * p_i;
%         
%         A_Rec = A_temp;
%         B_Rec = B_temp;
% 
%         X_test = [X_test x];
%         X_recon_test = [X_recon_test x_rec];
% 
%         disp(i / N * 100)
%     end
% 
%     figure
%     subplot(311)
%     plot(X_test(1, :))
%     legend('Original')
% 
%     subplot(312)
%     plot(X_test(2, :))
%     legend('Original')
% 
%     hold on 
%  
%     subplot(313)
% 
%     mesh(Tspan, 1 : Nd, P_i)
%     xlabel('Time(sec)')
%     ylabel('i-th vertex')
%     zlabel('p_i')
% 
%     figure
%     plot(U(:))
%     legend('Control input')
% 
% end
% 
% 















