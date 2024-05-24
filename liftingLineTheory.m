function [CL, CL1, y_s, CL_overall] = liftingLineTheory(N, alpha, alpha_0, span, C_r, C_L_alpha, t, AR)
    theta = pi/(2*N):pi/(2*N):pi/2;
    z = (span/2)*cos(theta);
    c = C_r * (1 - (1-t)*cos(theta)); % Mean Aerodynamics Chord at each segment (m)
    mu = c * C_L_alpha / (4 * span);
    LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
    B = zeros(N, N);
    for i=1:N
        for j=1:N
            B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
        end
    end
    A = B \ transpose(LHS);
    sum1 = zeros(1, N);
    sum2 = zeros(1, N);
    for i = 1:N
        for j=1:N
            sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
            sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
        end
    end
    CL = 4*span*sum2 ./ c;
    CL1 = [0 CL(1:N)];
    y_s = [span/2 z(1:N)];
    CL_overall = pi * AR * A(1);
end