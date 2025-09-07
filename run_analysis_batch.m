function run_analysis_batch(batch_def, batch_folder_name)

    %% 1. 파라미터 및 시뮬레이션 설정
    N = 5;
    bifurcation_param = batch_def.bifurcation_param;
    bifurcation_range = batch_def.bifurcation_range;
    fixed_params = batch_def.fixed_params;
    
    dt = 0.01;
    t_span_transient = 0:dt:20; % 과도 상태 제거용
    t_span_main = 0:dt:20;      % 분기점 수집용
    t_span_long = 0:dt:1000;     % 시간 구간 플로팅용

    rng(2024);
    x0 = rand(N,1);
    y0 = rand(N,1);
    ext_force_amp = 5.0;
    ext_force_freq = pi;

    %% 2. 분기 다이어그램 데이터 병렬 계산
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool; end
    
    num_steps = length(bifurcation_range);
    bifurcation_data = cell(num_steps, 1);
    
    fprintf('"%s" 파라미터에 대한 %d개 스텝의 병렬 계산을 시작합니다...\n', bifurcation_param, num_steps);
    parfor i = 1:num_steps
        current_val = bifurcation_range(i);
        
        p_vec = zeros(1, 4);
        fixed_names = fieldnames(fixed_params);
        for k = 1:length(fixed_names)
            param_name = fixed_names{k};
            switch param_name
                case 'a', p_vec(1) = fixed_params.a;
                case 'b', p_vec(2) = fixed_params.b;
                case 'c', p_vec(3) = fixed_params.c;
                case 'd', p_vec(4) = fixed_params.d;
            end
        end
        switch bifurcation_param
            case 'a', p_vec(1) = current_val;
            case 'b', p_vec(2) = current_val;
            case 'c', p_vec(3) = current_val;
            case 'd', p_vec(4) = current_val;
        end
        
        bifurcation_data{i} = calculate_bifurcation_point(p_vec, N, dt, t_span_transient, t_span_main, x0, y0, ext_force_amp, ext_force_freq);
    end
    fprintf('병렬 계산 완료.\n');

    %% 3. 분기 다이어그램 생성 및 저장
    bifurcation_folder = fullfile(batch_folder_name, 'Bifurcation_Diagrams');
    if ~exist(bifurcation_folder, 'dir'), mkdir(bifurcation_folder); end
    
    for nodeIdx = 1:N
        fig = figure('Visible','off');
        hold on;
        
        plot_bifurcation_scenario(bifurcation_data, bifurcation_range, nodeIdx, 'right', 'r.');
        plot_bifurcation_scenario(bifurcation_data, bifurcation_range, nodeIdx, 'mode', 'k.');
        plot_bifurcation_scenario(bifurcation_data, bifurcation_range, nodeIdx, 'left', 'b.');
        
        hold off;
        title(sprintf('Bifurcation Diagram (Node %d)', nodeIdx));
        xlabel(sprintf('Parameter %s', bifurcation_param));
        ylabel('Local Maxima of x');
        legend({'Right', 'Mode', 'Left'}, 'Location', 'northwest');
        grid on;
        
        file_path = fullfile(bifurcation_folder, sprintf('Bifurcation_Node%d.png', nodeIdx));
        exportgraphics(fig, file_path, 'Resolution', 150);
        close(fig);
    end
    fprintf('분기 다이어그램이 "%s" 폴더에 저장되었습니다.\n', bifurcation_folder);

    %% 4. 시간 구간별 그림 생성을 위한 최종 시뮬레이션
    fprintf('시간 구간별 그림 생성을 위한 최종 시뮬레이션을 실행합니다...\n');
    
    % --- 수정된 부분: parfor 루프 밖에서 final_p_vec 재생성 ---
    final_p_vec = zeros(1, 4);
    fixed_names = fieldnames(fixed_params);
    for k = 1:length(fixed_names)
        param_name = fixed_names{k};
        switch param_name
            case 'a', final_p_vec(1) = fixed_params.a;
            case 'b', final_p_vec(2) = fixed_params.b;
            case 'c', final_p_vec(3) = fixed_params.c;
            case 'd', final_p_vec(4) = fixed_params.d;
        end
    end
    % 분기 다이어그램의 마지막 파라미터 값을 사용
    last_bifurcation_value = bifurcation_range(end);
    switch bifurcation_param
        case 'a', final_p_vec(1) = last_bifurcation_value;
        case 'b', final_p_vec(2) = last_bifurcation_value;
        case 'c', final_p_vec(3) = last_bifurcation_value;
        case 'd', final_p_vec(4) = last_bifurcation_value;
    end
    % -----------------------------------------------------------

    numSteps_long = numel(t_span_long);
    [X_final, Y_final] = runSingleSimulation(dt, numSteps_long, x0, y0, final_p_vec, ext_force_amp, ext_force_freq, N);
    
    %% 5. 시간 구간별 그림 생성 및 저장
    interval_folder = fullfile(batch_folder_name, 'Interval_Phase_Portraits');
    plot_intervals(X_final, Y_final, final_p_vec, N, dt, interval_folder);
    fprintf('시간 구간별 그림이 "%s" 폴더에 저장되었습니다.\n', interval_folder);
end

function plot_bifurcation_scenario(data, range, nodeIdx, scenario, style)
    % 분기 다이어그램의 한 시나리오를 그리는 헬퍼 함수
    num_steps = length(range);
    all_peaks_x = [];
    all_peaks_y = [];
    
    for i = 1:num_steps
        if isfield(data{i}, scenario)
            peaks = data{i}.(scenario){nodeIdx};
            if ~isempty(peaks)
                all_peaks_x = [all_peaks_x, repmat(range(i), 1, length(peaks))];
                all_peaks_y = [all_peaks_y, peaks(:)'];
            end
        end
    end
    
    plot(all_peaks_x, all_peaks_y, style, 'MarkerSize', 2);
end
