clear;clc;close all

%% System parameters
A = [0.8,0.5;0.5,0.6];
B_fun = @(x) [0.5;x];
Q = eye(2); R = 0.1;

%% Constraint sets
lb_x = [-10,-5,0.5]; ub_x = [5,2,3];
Xc_vertex = [ub_x(1), ub_x(2);...
    ub_x(1), lb_x(2);...
    lb_x(1), lb_x(2);...
    lb_x(1), ub_x(2)];
Uc_vertex = [1; -1];
Xc = Polyhedron(Xc_vertex);
Uc = Polyhedron(Uc_vertex);

%% Generate samples
% Generate random samples inside the polyhedron
num_samples = 300;
[X,X_fea_family] = SampleGeneratorCCD(num_samples,lb_x,ub_x,A,B_fun,Uc_vertex);

scatter(X(:,1),X(:,2))

leg = legend({'$X_{fea}$', '$X_{fea}\ominus Z$', '$X_f$', '$X_f\ominus Z$','Samples'}, 'Location', 'southeast');
set(leg, 'Interpreter', 'latex','FontSize', 14);

%% Projection
% for any state x in Xfea, we project the network output uk onto a polytope
% (the set is to make sure the next state is still in the feasible set)
% Given u0 evaluated from the policy, the projecting is to yield u_p that
% satisfies the input constraints and the resulting x_next satisfies 
% maximal control invariant set
QP_setting.options = optimoptions('quadprog','Display','off');
QP_setting.f = @(u0) -2*u0;
QP_setting.H = 2;
QP_setting.lb_u = Uc_vertex(2);
QP_setting.ub_u = Uc_vertex(1);
for i = 1:num_samples
    Cc = X_fea_family{i}.A; dc = X_fea_family{i}.b;
    QP_setting.b_qp = @(x) dc-Cc*A*x;
    B = B_fun(X(i,end));
    QP_setting.A_qp = Cc*B;
    % Please provide u0
    x = X(i,1:2)';
    [K,S] = dlqr(A,B,Q,R);
    K = -K;
    u0 = K*x;
    % Solved by QP
    u_p_tmp = quadprog(QP_setting.H,QP_setting.f(u0),QP_setting.A_qp,QP_setting.b_qp(x),...
        [],[],QP_setting.lb_u,QP_setting.ub_u,u0,QP_setting.options);
    u_p(i) = u_p_tmp;
end
u_p = u_p(:);

figure;
scatter3(X(:,1),X(:,2),u_p',20,X(:,3),'filled');
colorbar