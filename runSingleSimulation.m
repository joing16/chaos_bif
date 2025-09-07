% 파일 이름: runSingleSimulation.m
function [X_full, Y_full] = runSingleSimulation(dt, numSteps, x0, y0, p_vec, ext_force_amp, ext_force_freq, N)
    phase_shift = 0; % 'mode' 시나리오 고정

    % 파라미터 행렬 생성
    [a, b, c, d] = deal(p_vec(1), p_vec(2), p_vec(3), p_vec(4));
    rng(2024);
    epsilon = 0.1;
    P = eye(N) + (rand(N) * 2 * epsilon - epsilon); P(logical(eye(N))) = 1;
    Q = eye(N) + (rand(N) * 2 * epsilon - epsilon); Q(logical(eye(N))) = 1;
    W = eye(N) + (rand(N) * 2 * epsilon - epsilon); W(logical(eye(N))) = 1;
    Z = eye(N) + (rand(N) * 2 * epsilon - epsilon); Z(logical(eye(N))) = 1;
    P = normalizeColumns(P); Q = normalizeColumns(Q);
    W = normalizeColumns(W); Z = normalizeColumns(Z);
    alpha_x = a * Q; beta_x = b * W;
    alpha_y = c * P; beta_y = d * Z;

    % CUDA 커널 호출
    [X_full, Y_full] = simulateModel2_full(dt, numSteps, x0, y0, alpha_x, beta_x, alpha_y, beta_y, ext_force_amp, ext_force_freq, N, phase_shift);
end
function Mnorm = normalizeColumns(M)
    for j = 1:size(M,2)
        colSum = sum(abs(M(:,j)));
        if colSum == 0, M(:,j) = 1/size(M,1); else, M(:,j) = M(:,j) / colSum; end
    end
    Mnorm = M;
end