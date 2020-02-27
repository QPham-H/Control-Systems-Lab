function funcWnZeta = parmSolve2ndOrderTF(x,N,T)
    % x = (zeta, wn)
    funcWnZeta = [x(1)*x(2) + log(N)/T;
                  sqrt(1-x(1)^2)*x(2) - 2*pi/T];
end