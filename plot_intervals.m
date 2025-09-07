%파일 이름: plot_intervals.m
function plot_intervals(X_full, Y_full, p_vec, N, dt, output_folder)
    intervals = {[21, 200], [221, 400], [421, 600], [621, 800], [821, 1000]};
    colors = {'b', 'g', 'r', 'c', 'm'};
    legend_labels = {'21-200s', '221-400s', '421-600s', '621-800s', '821-1000s'};
    
    % 제목 계산용 파라미터 행렬 생성
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

    for nodeIdx = 1:N
        node_output_folder = fullfile(output_folder, sprintf('Node_%d', nodeIdx));
        if ~exist(node_output_folder, 'dir'), mkdir(node_output_folder); end

        sum_ax = sum(alpha_x(:, nodeIdx)); sum_bx = sum(beta_x(:, nodeIdx));
        sum_ay = sum(alpha_y(:, nodeIdx)); sum_by = sum(beta_y(:, nodeIdx));
        paramTitle = sprintf('α_x=%.2f, β_x=%.2f, α_y=%.2f, β_y=%.2f', sum_ax, sum_bx, sum_ay, sum_by);

        % 구간별 개별 그림
        for i = 1:length(intervals)
            startTime = intervals{i}(1); endTime = intervals{i}(2);
            startIndex = floor(startTime / dt) + 1; endIndex = floor(endTime / dt) + 1;

            fig = figure('Visible','off');
            x_data = X_full(nodeIdx, startIndex:endIndex); y_data = Y_full(nodeIdx, startIndex:endIndex);
            finite_idx = isfinite(x_data) & isfinite(y_data);
            plot(x_data(finite_idx), y_data(finite_idx), 'Color', colors{i});
            title({sprintf('Phase Portrait (%s)', legend_labels{i}), paramTitle});
            xlabel(sprintf('x_%d', nodeIdx)); ylabel(sprintf('y_%d', nodeIdx)); grid on;
            
            fileName = fullfile(node_output_folder, sprintf('Phase_%s.png', legend_labels{i}));
            exportgraphics(fig, fileName, 'Resolution', 150); close(fig);
        end
        
        % 통합 그림
        fig_total = figure('Visible','off');
        hold on;
        for i = 1:length(intervals)
            startTime = intervals{i}(1); endTime = intervals{i}(2);
            startIndex = floor(startTime / dt) + 1; endIndex = floor(endTime / dt) + 1;
            x_data = X_full(nodeIdx, startIndex:endIndex); y_data = Y_full(nodeIdx, startIndex:endIndex);
            finite_idx = isfinite(x_data) & isfinite(y_data);
            plot(x_data(finite_idx), y_data(finite_idx), 'Color', colors{i});
        end
        hold off;
        title({'Total Phase Portrait', paramTitle});
        xlabel(sprintf('x_%d', nodeIdx)); ylabel(sprintf('y_%d', nodeIdx)); grid on;
        legend(legend_labels, 'Location', 'southeast');
        
        fileName = fullfile(node_output_folder, 'Phase_Total.png');
        exportgraphics(fig_total, fileName, 'Resolution', 150); close(fig_total);
    end
end

function Mnorm = normalizeColumns(M)
    for j = 1:size(M,2)
        colSum = sum(abs(M(:,j)));
        if colSum == 0, M(:,j) = 1/size(M,1); else, M(:,j) = M(:,j) / colSum; end
    end
    Mnorm = M;
end
