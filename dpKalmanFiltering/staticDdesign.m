clear;
clc;

%% If need random numbers with repeatable experiments
%rng('default');
%rng(12345);

%% Problem data: A, W, C, V for agent i
Ts = 1;
A = [1 Ts;0 1];
%A = [0 1;0 1];

n = size(A,1);

% B = [Ts^2/2;Ts];
% %sig = 1;
% %W = B*sig*B';
% %W = eye(2);
% W = 0.5*randn(2); W = W*W';
W = 1e-1 * eye(n);

C = eye(n);
m = size(C,1);

Winv = inv(W);
V = eye(m);

%% Solve the SDP
alpha = 0.8;
pibound = inv(alpha^2 * eye(m) + V);

cvx_begin sdp
    variable S(n,n) semidefinite;  % Sigma matrix, psd
    variable R(n,n) semidefinite;  % Omega matrix, inverse of Sigma, psd
    variable P(m,m) semidefinite;  % Pi matrix, psd
    [S, eye(n); eye(n), R] >= 0;
    
    %[R-C'*P*C, eye(n); eye(n), A*S*A'+W] >= 0;
    [C'*P*C-R+Winv, Winv*A; A'*Winv, R+A'*Winv*A] >= 0;
    
    P <= pibound;
    
    minimize(trace(S));
cvx_end

disp('Difference between inv(S) and R')
E1 = inv(S)-R
norm(E1)

disp('Difference between inv(S) and inv(A S A^T + W)+P')
E2 = inv(S)-inv(A*S*A'+W)-P
norm(E2)

disp('Matrix to factorize')
M1 = alpha^2 * (inv(V-V*P*V)-inv(V))

%M = [S,eye(n);eye(n),R];
%eig(M)

%S
%inv(R)