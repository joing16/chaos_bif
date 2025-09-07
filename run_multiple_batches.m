function run_multiple_batches()
    clear; clc;
    
    fprintf('분기 다이어그램 및 시간 구간 분석 자동 실행을 시작합니다...\n\n');
    
    % =====================================================================
    % 여기에 실행할 모든 분석 회차를 미리 정의합니다.
    % =====================================================================
    batch_definitions = {
        % --- 회차 1 ---
        struct('bifurcation_param', 'a', ...
               'bifurcation_range', linspace(-30, -0.1, 500), ...
               'fixed_params', struct('b', -1.0, 'c', 0.86, 'd', 0.5)),
       %  % B바
       % 
       %  struct('bifurcation_param', 'b', ...
       %         'bifurcation_range', linspace(-0.005, -3, 20), ...
       %         'fixed_params', struct('a', -3.0, 'c', 0.86, 'd', 0.5)),
       % 
       % %C바               
       %  struct('bifurcation_param', 'c', ...
       %         'bifurcation_range', linspace(0.01, 3, 20), ...
       %         'fixed_params', struct('a', -3.0, 'b', -1.0, 'd', 0.5)),
       % % D바
       %  struct('bifurcation_param', 'd', ...
       %         'bifurcation_range', linspace(0.0001, 1, 20), ...
       %         'fixed_params', struct('a', -3.0, 'b', -1.0, 'c', 0.86)),
       % 
        % 여기에 필요한 만큼 회차(struct)를 계속 추가할 수 있습니다.
    };
    % =====================================================================

    num_batches = length(batch_definitions);
    
    for i = 1:num_batches
        fprintf('====================================================\n');
        fprintf('>>> 회차 %d / %d 실행 시작 (시간: %s)\n', i, num_batches, datestr(now,'HH:MM:SS'));
        fprintf('====================================================\n');
        
        current_batch = batch_definitions{i};
        
        % 현재 회차의 파라미터 정보로 고유 폴더 이름 생성
        fp = current_batch.fixed_params;
        param_names = fieldnames(fp);
        folder_str = '';
        for k = 1:length(param_names)
            folder_str = [folder_str, sprintf('%s%.2f_', param_names{k}, fp.(param_names{k}))];
        end
        batch_folder_name = sprintf('Batch%d_Bifurcation(%s)_%s', i, current_batch.bifurcation_param, folder_str);
        batch_folder_name = strrep(batch_folder_name, '.', '_');
        
        % 분석 실행 함수 호출
        run_analysis_batch(current_batch, batch_folder_name);
                               
        fprintf('\n>>> 회차 %d / %d 실행 완료 (시간: %s)\n\n', i, num_batches, datestr(now,'HH:MM:SS'));
    end
    
    fprintf('모든 회차의 분석이 완료되었습니다.\n');
end
