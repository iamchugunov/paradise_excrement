function [X, dop, nev, flag] = NavSolver_D(y, posts, X0, h)

    epsilon = 0.001;
    max_iter = 7;
    
    N = size(posts,2);
    Y = zeros(N, 1);
    H = zeros(N, 2);
    X = X0;
    
    iter = 0;
    
while 1
    
    iter = iter + 1;
    
    for i = 1:N
        d = sqrt((X(1,1) - posts(1,i))^2 + (X(2,1) - posts(2,i))^2 + (h - posts(3,i))^2);
        Y(i, 1) = d;
        H(i, 1) = (X(1,1) - posts(1,i))/d;
        H(i, 2) = (X(2,1) - posts(2,i))/d;        
    end
    
    X_prev = X;
    X = X + (H'*H)^(-1)*(H')*(y-Y);
    
    if (norm(X - X_prev) < epsilon) || (iter > max_iter)
        invHH = inv(H'*H);
        DOPx = abs(invHH(1,1));
        DOPy = abs(invHH(2,2));
        dop = sqrt([DOPx DOPy])';
        break
    end
end

if norm(X - X_prev) < 1e5
    flag = 1;
    nev = norm(y - Y);
    if nev > 0.4
        flag = 0;
    end
else
    flag = 0;
    nev = 0;
end

end

