%% Task 1: State-space model
syms sgma omega x1 x2 x3 x4 u1 u2

% 1.1 Provide the following expressions with symbolic variables
g=omega^2*sgma^3;
x1_dot=x2
x2_dot=u1 - g/(x1 + sgma)^2 + (x1 + sgma)*(x4/sgma + omega)^2
x3_dot=x4
x4_dot=sgma/(x1 + sgma)*(u2 - 2*x2*(x4/sgma + omega))

x_dot = [x1_dot;x2_dot;x3_dot;x4_dot];
% State space matrices in symbolic expressions
A_sym=jacobian(x_dot, [x1, x2, x3, x4])
B_sym=jacobian(x_dot, [u1, u2])


% 1.2 Let sgma and omega all equal to 1. Provide the numeric matrices of the state-space system, in double-precision floating-point format
A=subs(A_sym, [sgma omega x1 x2 x3 x4 u1 u2], [1, 1, 0, 0, 0, 0, 0, 0])
B=subs(B_sym, [sgma omega x1 x2 x3 x4 u1 u2], [1, 1, 0, 0, 0, 0, 0, 0])
C=[1 0 0 0; 0 0 1 0]
D=zeros(2)

%% Task 2: Reachability
% 2.1 
W_r=[B, A*B, A^2*B, A^3*B]
length(A) - rank(W_r)
Sys_reachable=1 % equal to 1 if system reachable, equal to 0 if not.

% 2.2
W_r_1=[B(:, 2), A*B(:, 2), A^2*B(:, 2), A^3*B(:, 2)]  %when u1=0
length(A) - rank(W_r_1)
Sys_reachable_1=1 % equal to 1 if system reachable, equal to 0 if not.

% 2.3
W_r_2=[B(:, 1), A*B(:, 1), A^2*B(:, 1), A^3*B(:, 1)] %when u2=0
length(A) - rank(W_r_2)
Sys_reachable_2=0 % equal to 1 if system reachable, equal to 0 if not.