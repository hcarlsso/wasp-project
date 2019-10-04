w = randn(3,1);
w_dot = randn(3,1);
s = randn(3,1);
Na = 3;
Ng = 1;
r = randn(3, Na);

ya = zeros(3, Na);
yg = zeros(3, Ng);

sig = 0.01;
for k = 1:Na
    ya(:,k) = s + cross(w, cross(w, r(:,k))) + cross(w_dot, r(:,k)) + sig*randn(3,1);
end

for k = 1:Ng
    yg(:,k) = w + sig*randn(3,1);
end

addpath('solvers/')

[s_hat, w_hat, w_dot_hat] = solveImuArray(ya, yg, r,1.0,1.0);

disp('s')
s_hat - s

disp('w')
w_hat-w

disp('w_dot')
w_dot - w_dot_hat