function [H_traj, K_traj, A_traj, HEEL_traj, TOE_traj] = compute_trajectories(q_cmd, params)
    % 主函数：批量2D正运动学轨迹计算（对应论文 Results 足部轨迹图）
    % 输入：q_cmd (N x 3) 关节角度序列，params 结构体
    % 输出：各关节点轨迹矩阵 (N x 2)
    
    % 输入检查
    assert(size(q_cmd,2)==3, 'q_cmd 必须为 N×3 矩阵');
    
    N = size(q_cmd, 1);
    model = params.model;
    
    % 预分配数组（提升性能）
    H_traj    = zeros(N, 2);
    K_traj    = zeros(N, 2);
    A_traj    = zeros(N, 2);
    HEEL_traj = zeros(N, 2);
    TOE_traj  = zeros(N, 2);
    
    for i = 1:N
        q1 = q_cmd(i,1); 
        q2 = q_cmd(i,2); 
        q3 = q_cmd(i,3);
        
        pts = forward_kinematics_2d(q1, q2, q3, model);
        
        H_traj(i,:)    = pts.H';
        K_traj(i,:)    = pts.K';
        A_traj(i,:)    = pts.A';
        HEEL_traj(i,:) = pts.HEEL';
        TOE_traj(i,:)  = pts.TOE';
    end
    
    % 可选：轻微平滑（使论文图表更连续美观）
    % TOE_traj = smoothdata(TOE_traj, 1, 'movmean', 5);
end

function pts = forward_kinematics_2d(q1, q2, q3, model)
    % 单时刻2D正运动学（L型脚建模）
    H = [0; 0];

    theta1 = q1;
    theta2 = q1 + q2;
    theta3 = q1 + q2 + q3;

    K = H + [model.L1*cos(theta1); model.L1*sin(theta1)];
    A = K + [model.L2*cos(theta2); model.L2*sin(theta2)];

    HEEL = A + [model.L_ankle_to_heel*cos(theta3); model.L_ankle_to_heel*sin(theta3)];
    TOE  = HEEL + [model.L_heel_to_toe*cos(theta3 + pi/2); model.L_heel_to_toe*sin(theta3 + pi/2)];

    pts.H    = H;
    pts.K    = K;
    pts.A    = A;
    pts.HEEL = HEEL;
    pts.TOE  = TOE;
end