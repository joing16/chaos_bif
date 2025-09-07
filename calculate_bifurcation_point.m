function data_point = calculate_bifurcation_point(p_vec, N, dt, t_span_transient, t_span_main, x0, y0, ext_force_amp, ext_force_freq)
    
    scenarios = {'right', 'mode', 'left'};
    data_point = struct();

    % 파라미터 행렬 생성
    [a, b, c, d] = deal(p_vec(1), p_vec(2), p_vec(3), p_vec(4));
    seed = sum(uint32(mat2str(p_vec))) + N;
    rng(seed);
    epsilon = 0.1;
    P = eye(N) + (rand(N) * 2 * epsilon - epsilon); P(logical(eye(N))) = 1;
    Q = eye(N) + (rand(N) * 2 * epsilon - epsilon); Q(logical(eye(N))) = 1;
    W = eye(N) + (rand(N) * 2 * epsilon - epsilon); W(logical(eye(N))) = 1;
    Z = eye(N) + (rand(N) * 2 * epsilon - epsilon); Z(logical(eye(N))) = 1;
    P = normalizeColumns(P); Q = normalizeColumns(Q);
    W = normalizeColumns(W); Z = normalizeColumns(Z);
    alpha_x = a * Q; beta_x = b * W;
    alpha_y = c * P; beta_y = d * Z;

    % 각 시나리오에 대해 계산
    for i = 1:length(scenarios)
        scenario = scenarios{i};
        
        switch scenario
            case 'right', phase_shift = 0.2;
            case 'mode',  phase_shift = 0;
            case 'left',  phase_shift = -0.2;
        end
        
        % 1. 과도 상태(transient) 제거를 위한 초기 시뮬레이션
        numSteps_transient = numel(t_span_transient);
        [X_trans, Y_trans] = simulateModel2_full(dt, numSteps_transient, x0, y0, alpha_x, beta_x, alpha_y, beta_y, ext_force_amp, ext_force_freq, N, phase_shift);
        
        % 초기 시뮬레이션의 마지막 상태를 본 시뮬레이션의 초기값으로 사용
        x0_main = X_trans(:, end);
        y0_main = Y_trans(:, end);
        
        % 2. 분기점 수집을 위한 본 시뮬레이션
        numSteps_main = numel(t_span_main);
        [X_main, ~] = simulateModel2_full(dt, numSteps_main, x0_main, y0_main, alpha_x, beta_x, alpha_y, beta_y, ext_force_amp, ext_force_freq, N, phase_shift);
        
        % 3. 각 노드의 극대값(peaks) 찾기
        node_peaks = cell(N, 1);
        for nodeIdx = 1:N
            x_series = X_main(nodeIdx, :);
            finite_idx = isfinite(x_series);
            
            % --- 수정된 부분: findpeaks 호출 전 데이터가 비어있는지 확인 ---
            valid_x_series = x_series(finite_idx);
            if ~isempty(valid_x_series)
                [peaks, ~] = findpeaks(valid_x_series);
            else
                peaks = []; % 유효한 데이터가 없으면 peaks를 빈 배열로 설정
            end
            % -----------------------------------------------------------
            
            node_peaks{nodeIdx} = peaks;
        end
        
        data_point.(scenario) = node_peaks;
    end
end

function Mnorm = normalizeColumns(M)
    for j = 1:size(M,2)
        colSum = sum(abs(M(:,j)));
        if colSum == 0, M(:,j) = 1/size(M,1); else, M(:,j) = M(:,j) / colSum; end
    end
    Mnorm = M;
end
