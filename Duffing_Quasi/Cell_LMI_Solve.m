function [g, u, Q_tilde] = Cell_LMI_Solve(Now_A, Now_B, A, B, xk, Nd)
    n = 4;
    m = 1;
    C = [1 0 0 0];
    Q = 10 * diag([1, 1, 0, 0]);
    R = 0.01;
    y_max = 10;
    u_max = 4;
    g = sdpvar(1, 1);
    u = sdpvar(m, 1);
    X = sdpvar(m, m);
    Q_tilde = sdpvar(n, n);
%     xk = zeros(n, 1);
    M = {};
    w = 0.17* ones(n, 1);
    epsilon1 = 1;
    epsilon2 = 1 / epsilon1;
    Input_C = {};
    Output_C = {};
    Y = sdpvar(m, n);
    M1 = [1 (Now_A*xk+Now_B*u)'*sqrt(1 + epsilon1) xk'*sqrt(Q) u'*sqrt(R) sqrt(1 + epsilon2) * w';
          (Now_A*xk+Now_B*u)*sqrt(1 + epsilon1) Q_tilde zeros(n, n) zeros(n, m) zeros(n, n);
          sqrt(Q)*xk zeros(n, n) (g)*eye(n) zeros(n, m) zeros(n, n);
          sqrt(R)*u zeros(m, n) zeros(m, n) (g )*eye(m) zeros(m, n);
          sqrt(1 + epsilon2) * w zeros(n, n) zeros(n, n) zeros(n, m)  Q_tilde];
%     M1 = [1 (Now_A*xk+Now_B*u)' xk'*sqrt(Q) u'*sqrt(R);
%           Now_A*xk+Now_B*u Q_tilde zeros(n, n) zeros(n, m);
%           sqrt(Q)*xk zeros(n, n) g*eye(n) zeros(n, m);
%           sqrt(R)*u zeros(m, n) zeros(m, n) g*eye(m)];
    for i = 1 : Nd
        M{i} = [Q_tilde Q_tilde*A{i}'+Y'*B{i}' Q_tilde*sqrt(Q) Y'*sqrt(R);
          A{i}*Q_tilde+B{i}*Y Q_tilde zeros(n, n) zeros(n, m);
          sqrt(Q)*Q_tilde zeros(n, n) g*eye(n) zeros(n, m);
          sqrt(R)*Y zeros(m, n) zeros(m, n) g*eye(m)];
        Input_C{i} = [X Y; Y' Q_tilde];
        Output_C{i} = [Q_tilde, (A{i} * Q_tilde + B{i} * Y)'*C';
                 C*(A{i} * Q_tilde + B{i} * Y) y_max^2];
    end

    Input_C1 = [u - u_max; 
               -u_max - u];
    Output_C1 = [C*(Now_A*xk + Now_B*u) - y_max;
                -y_max - C*(Now_A*xk + Now_B*u)];
    Constraints = [M1 >= 0, Input_C1 <= 0];
    for i = 1 : Nd
        Constraints = [Constraints, M{i} >= 0, Input_C{i} >= 0];
    end
    for i = 1 : m
        Constraints = [Constraints X(i, i) <= u_max^2];
    end

    obj = g;
    
    sol = solvesdp(Constraints,obj)
    if sol.problem ~= 0 
        error('infeasible')
    end
    X = double(X)
    g = double(g)
    u = double(u)
    Q_tilde = double(Q_tilde)
    
    
    
end





