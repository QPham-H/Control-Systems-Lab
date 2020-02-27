function M = MCal(zeta,wn,t)
    M = 1/sqrt(1-zeta^2)*exp(-zeta*wn*t)*sin(wn*sqrt(1-zeta^2)*t+acos(zeta));
end