w = randn(3,1)*20;
w_dot = randn(3,1)*100;
s = randn(3,1);
% Na = 3; r = randn(3, Na);
Ng = 0;

r = [0 1 0 0 
     0 0 1 0 
     0 0 0 1 ];
    
Na = size(r,2);

ya = zeros(3, Na);
yg = zeros(3, Ng);

sig = 0.0;
for k = 1:Na
    ya(:,k) = s + cross(w, cross(w, r(:,k))) + cross(w_dot, r(:,k)) + sig*randn(3,1);
end

for k = 1:Ng
    yg(:,k) = w + sig*randn(3,1);
end

addpath('solvers/')

[s_hat, w_hat, w_dot_hat, L] = solveImuArray(ya, yg, r, 1.0,1.0);
if Ng > 0
    w0 = yg(:,1);
else
    w0 = w ;
end
w0 = randn(3,1);
w_gs = solveImuArrayGs(ya,yg,r,1.0,1.0, w0);

disp('results w')
w
w_gs
w_hat(:,1)
w_hat
L

if true
    %% 
    wx_range = linspace(min(w_hat(1,:))-0.5,  max(w_hat(1,:))+0.5, 100);
    wy_range = linspace(min(w_hat(2,:))-0.5,  max(w_hat(2,:))+0.5, 100);
    wz_range = linspace(min(w_hat(3,:))-0.5,  max(w_hat(3,:))+0.5, 20);

    [Wx, Wy, Wz] = meshgrid(wx_range, wy_range, wz_range);


    cost = concentrated_cost(ya,yg, r, 1.0, 1.0, Wx, Wy, Wz);

    %% 
    figure(1); clf
    % scatter3(Wx(:), Wy(:), Wx(:), 1.0, cost(:))
    xslice = [];   
    yslice = [];
    zslice = wz_range;
    lvls = logspace(-3, log10(max(cost(:))), 30);
    contourslice(Wx, Wy, Wz, cost,xslice,yslice,zslice, lvls)

    hold on 
    plot3(w(1),w(2),w(3), 'rs', 'MarkerSize', 10)

    plot3(w_hat(1,1), w_hat(2,1), w_hat(3,1), 'rx', 'MarkerSize', 12)
    S = -1./(L-L(1) - 1.5e-1) + 10;
    scatter3(w_hat(1,2:end), w_hat(2,2:end),w_hat(3,2:end), S(2:end) ,"ro")
    plot3(w_gs(1,:), w_gs(2,:), w_gs(3,:), '-r')
    colorbar
    xlabel("x")
    ylabel("y")
    zlabel("z")

    view(3)
    grid on

    %%
end