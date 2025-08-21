function  [Alift, Blift, Clift, Ylift, Xlift, X, U] = Duffing_System_Generate (liftFun,Ubig, Xcurrent, Nsim, Ntraj)


%% *************************** Dynamics ***********************************
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 3*x(1,:).^2.*x(2,:) - 1*x(2,:) - u] );
% f_u =  @(t,x,u)([ x(2,:) ; 4*9.8*sin(x(1, :)) - u ] );
% f_u =  @(t,x,u)([ x(2,:) ; x(1, :) - x(2, :).^3 + u] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) - u]);
% f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) - u]);
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)([ x(2,:) ; x(1, :) - x(2, :).^3 + u] );
f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );

% f_u =  @(t,x,u)(-[ -2*x(2,:) ; alpha * x(1,:) + beta*x(1,:).^2.*x(2,:) - gamma*x(2,:) + u] );
% f_u =  @(t,x,u)([ x(2,:) ; -alpha*x(2, :) + beta * x(1, :) - gamma * x(1, :).^3 + u] );
% f_u =  @(t,x,u)([ x(2,:) ; 4*9.8*sin(x(1, :)) - u .* cos(x(1, :))] );
% f_u =  @(t,x,u)([ x(2,:) ;x(1, :).^2+ 0.15 .* u(1, :).^3  + 0.1* (1 + x(2,:).^2).*u(1, :) + sin(0.1 .* u(1, :))] );

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

% basisFunction = 'rbf';
% % RBF centers
% Nrbf = 5;
% cent = rand(n,Nrbf)*2 - 1;
% rbf_type = 'gauss'; 
% % Lifting mapping - RBFs + the state itself
% % liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
% % liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] - [zeros(2, 1);rbf(zeros(2, 1),cent,rbf_type)]);
% 
% % liftFun = @(xx)( [xx; Encoder_VDP(xx)] - [zeros(2, 1); Encoder_VDP(zeros(2, 1))] );
% % Nlift = Nrbf + n;
% % liftFun = @(xx)( [xx]);
% liftFun = @(xx)( [xx; xx(1).*xx(2); xx(1)^2.*xx(2); xx(2)^2.*xx(1)]);
% liftFun = @(xx)( [xx; xx(1).*xx(2); xx(1)^2.*xx(2); xx(2)^2.*xx(1); xx(2)^2.*xx(1).^2]);

Nlift =size( liftFun(zeros(2, 1)), 1);


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
% if(exist('Nsim','var'))
%     Nsim = 200;  % 如果未出现该变量，则对其进行赋值
% end

% if(exist('Ntraj','var'))
%     Ntraj = 1000;  % 如果未出现该变量，则对其进行赋值
% end


% Random forcing
% Ubig = 2*rand([Nsim Ntraj]) - 1;
% if(exist('U1','var'))
%     Ubig = U1;  % 如果未出现该变量，则对其进行赋值
% end

% Random initial conditions
% if(exist('Xcurrent','var'))
%     Xcurrent = (rand(n,Ntraj)*2 - 1); % 如果未出现该变量，则对其进行赋值
% end

X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift = [];
Ylift = [];
for i = 1 : size(X, 2)
    Xlift = [Xlift liftFun(X(:, i))];
    Ylift = [Ylift liftFun(Y(:, i))];
end
% Xlift = liftFun(X);
% Ylift = liftFun(Y);
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic
W = [Ylift ; X];
V = [Xlift; U];
VVt = V*V';
WVt = W*V';
M = WVt * pinv(VVt); % Matrix [A B; C 0]
Alift = M(1:Nlift,1:Nlift);
Blift = M(1:Nlift,Nlift+1:end);
Clift = M(Nlift+1:end,1:Nlift);

fprintf('Regression done, time = %1.2f s \n', toc);
