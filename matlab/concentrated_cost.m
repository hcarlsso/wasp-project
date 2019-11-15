function cost = concentrated_cost(ya,yg,r,sa,sg, Wx, Wy, Wz)

    % Sanitate input data.
    ya = ya(:);
    yg = yg(:);
    y = [ya; yg];

    Na = length(ya)/3;
    Ng = length(yg)/3;

    assert(Na >= 3,'There has to be at least three accelerometers.');
    assert(isequal(size(r),[3 Na]),'The size of r is inconsitent with ya. Expected a size of [3 %d].',Na);

    % Create problem matrices.
    Ha = [-skewSymmetric(r) repmat(eye(3),Na,1)];

    Qai = (1/sa)*eye(3*Na);
    Qgi = (1/sg)*eye(3*Ng);

    HQ = (Ha'*Qai*Ha)\Ha'*Qai;
    Pa = Qai-Qai*Ha*HQ;
    Pg = Qgi;
    
    P = blkdiag(Pa,Pg);

    E = zeros(9,9);
    E([5 9 11 13 21 25 28 36 42 44 46 50]) = [-1 -1 1 1 1 1 -1 -1 1 1 -1 -1];
    Wa = kron(r',eye(3))*E;
    Wg = [zeros(3*Ng,6) repmat(eye(3),Ng,1)];
    W = [Wa; Wg]; % h(w) = W*m
    

    
    cost = zeros(size(Wx));

    for i = 1:size(cost,1)
        for j = 1:size(cost,2)
            for k = 1:size(cost,3)
                wx = Wx(i,j,k);
                wy = Wy(i,j,k);
                wz = Wz(i,j,k);
                
                m = [wx^2; wx*wy; wx*wz; wy^2; wy*wz; wz^2; wx; wy; wz];
                e = y-W*m;

                cost(i,j,k) = e'*P*e;
                
            end
        end
    end
end