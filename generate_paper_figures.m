%% generate_paper_figures.m — 四控制器论文 Fig.5～Fig.11 + supervisory diagnostic

function generate_paper_figures(t_all, q_ref, dq_ref, ddq_ref, ...
                                q_fixed, dq_fixed, ddq_fixed, ...
                                q_zhang, dq_zhang, ddq_zhang, ...
                                q_cao, dq_cao, ddq_cao, ...
                                q_adapt, dq_adapt, ddq_adapt, ...
                                Ba_zhang, Ka_zhang, Ga_zhang, ...
                                Ba_cao, Ka_cao, Ga_cao, ...
                                Ba_adp, Ka_adp, Ga_adp, ...
                                metrics_fixed, metrics_zhang, metrics_cao, metrics_adaptive, ...
                                params, sim_timestamp, monte_carlo_path)

    %% ---------- 0. 路径 ----------
    results_root = 'results';
    sim_folder   = fullfile(results_root, sim_timestamp);
    if ~exist(sim_folder, 'dir')
        error('仿真结果文件夹 %s 不存在', sim_folder);
    end

    fig_dir = fullfile(sim_folder, 'paper_figures');
    dat_dir = fullfile(sim_folder, 'figure_data');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    if ~exist(dat_dir, 'dir'),  mkdir(dat_dir); end

    
    if nargin < 32 || isempty(monte_carlo_path)
        mc_glob = fullfile(results_root, 'MC_*');
        mc_dirs = dir(mc_glob);
        if ~isempty(mc_dirs)
            [~, idx] = max([mc_dirs.datenum]);
            monte_carlo_path = fullfile(mc_dirs(idx).folder, mc_dirs(idx).name);
        else
            monte_carlo_path = '';
        end
    end

    %% ---------- 1. 全局样式 ----------
    % 论文图保留整图大标题，同时保留每个分图上方的小标题。
    % 布局采用显式坐标轴位置，为大标题、分图标题和每幅图的 x 轴刻度数字分别预留空间。
    SHOW_MAIN_TITLE = true;
    JPEG_DPI = 1000;  % journal-ready raster export requested for submission
    RASTER_PAD_PT = 18; % explicit white border around JPEG (18 pt = 0.25 in)

    set(groot, 'DefaultAxesFontSize', 11);
    set(groot, 'DefaultTextFontSize', 12);
    set(groot, 'DefaultLineLineWidth', 1.6);
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1);
    set(groot, 'DefaultAxesTitleFontWeight', 'normal');
    % 所有论文分图保留主网格线。设置为 root 默认值可防止 plot/bar/boxplot
    % 在首次绘制并重置坐标轴属性时把网格恢复为关闭状态。
    set(groot, 'DefaultAxesXGrid', 'on');
    set(groot, 'DefaultAxesYGrid', 'on');

    color_ref   = [0.2 0.2 0.2];        % 参考轨迹
    color_fixed = [0.0 0.4470 0.7410];  % Fixed
    color_zhang = [0.4660 0.6740 0.1880];  % Zhang-2018 VA
    color_cao   = [0.4940 0.1840 0.5560];  % Cao-2024-derived VA
    color_adapt = [0.8500 0.3250 0.0980];  % Proposed IE-AAC
    joints_label = {'Hip','Knee','Ankle'};
    N = length(t_all);
    dt = params.dt;

    % 论文图的事件窗口与真实扰动参数自动同步
    [~, event_timing] = event_model(t_all, params);

    %% ---------- 2. 辅助函数 ----------
    function export_fig(fig, name)
        % 重要：PDF 白边必须由 PDF page size 本身提供。
        % exportgraphics(fig,...) 会按图形内容/figure 边界输出，
        % 不会可靠继承 PaperPosition / PaperSize 中人为增加的纸张白边。
        % 因此论文 PDF 统一使用 print -dpdf -painters 导出。
        drawnow;
        set(fig, 'Renderer', 'painters', 'Color', 'w', 'InvertHardcopy', 'off');

        % -------- PDF 页面四周真实留白 --------
        % 四周使用 18 pt (= 0.25 in) 的明确白边，避免插入 Word/稿件后
        % 图中标题、坐标文字与正文在视觉上紧贴。
        pdf_pad_pt = 18;

        old_units = get(fig, 'Units');
        set(fig, 'Units', 'points');
        fig_pos_pt = get(fig, 'Position');
        fig_w_pt = max(fig_pos_pt(3), 1);
        fig_h_pt = max(fig_pos_pt(4), 1);

        set(fig, 'PaperUnits', 'points', ...
                 'PaperPositionMode', 'manual', ...
                 'PaperPosition', [pdf_pad_pt, pdf_pad_pt, fig_w_pt, fig_h_pt], ...
                 'PaperSize', [fig_w_pt + 2*pdf_pad_pt, ...
                               fig_h_pt + 2*pdf_pad_pt]);

        pdf_filename = fullfile(fig_dir, [name, '.pdf']);
        jpg_filename = fullfile(fig_dir, [name, '.jpg']);

        % -------- 1000 dpi JPEG + explicit outer white border --------
        % MATLAB R2025a+ supports exportgraphics(...,'Padding',...).
        % The default exportgraphics padding is tight, which is exactly what
        % caused the JPEG boundary to sit against titles/labels.  Request an
        % explicit 18-pt white border.  Older releases fall back to print;
        % print preserves the figure/page whitespace configured above.
        try
            exportgraphics(fig, jpg_filename, ...
                           'Resolution', JPEG_DPI, ...
                           'BackgroundColor', 'white', ...
                           'Units', 'points', ...
                           'Padding', RASTER_PAD_PT);
        catch
            try
                print(fig, jpg_filename, '-djpeg100', sprintf('-r%d', JPEG_DPI));
            catch
                try
                    print(fig, jpg_filename, '-djpeg', sprintf('-r%d', JPEG_DPI));
                catch MEjpg
                    set(fig, 'Units', old_units);
                    error('generate_paper_figures:JPEGExportFailed', ...
                          '1000 dpi JPEG 导出失败：%s', MEjpg.message);
                end
            end
        end

        try
            % print 会严格使用上面的自定义 PaperSize / PaperPosition，
            % 因而四边都能真正留下 pdf_pad_pt 的页面白边。
            print(fig, pdf_filename, '-dpdf', '-painters', '-r300');
        catch ME
            % 不再静默回退到 exportgraphics，因为它会重新产生“贴边”问题。
            set(fig, 'Units', old_units);
            error('generate_paper_figures:PDFExportFailed', ...
                  ['PDF 导出失败：%s\n' ...
                   '为保证论文 PDF 四周白边，本版本不回退到 exportgraphics。'], ...
                  ME.message);
        end

        set(fig, 'Units', old_units);
        fprintf('   ✅ 图 %s.pdf + %s.jpg 保存（JPEG %d dpi；JPEG/PDF 四周约 %.0f pt 白边）\n', ...
                name, name, JPEG_DPI, RASTER_PAD_PT);
    end

    function [fig, layout] = create_paper_layout(fig_name, fig_pos, nrow, ncol, main_title)
        % 使用显式归一化坐标布置各分图，而不再依赖 tiledlayout 自动压缩。
        % 这样可以同时满足：
        %   1) 大标题保留；
        %   2) 分图标题位于每个坐标轴正上方；
        %   3) 3x1 图中每一幅都保留 x 轴刻度数字；
        %   4) 不同长度的 y 轴刻度不会改变各分图绘图区的水平对齐。
        fig = figure('Name', fig_name, 'Position', fig_pos, 'Color', 'w');

        % 针对论文中三种主要版式分别预留足够的上下/左右空白。
        if nrow == 3 && ncol == 1
            left_margin   = 0.115;
            right_margin  = 0.045;
            bottom_margin = 0.075;
            top_margin    = 0.118;   % 收紧大标题与第一幅分图标题之间的距离
            row_gap       = 0.105;   % 上一幅 x 刻度数字 + 下一幅分图标题
            col_gap       = 0.0;
        elseif nrow == 1 && ncol == 3
            left_margin   = 0.070;
            right_margin  = 0.035;
            bottom_margin = 0.155;
            top_margin    = 0.118;   % 收紧大标题与三幅分图标题之间的距离
            row_gap       = 0.0;
            col_gap       = 0.080;
        elseif nrow == 3 && ncol == 3
            left_margin   = 0.070;
            right_margin  = 0.035;
            bottom_margin = 0.075;
            top_margin    = 0.118;   % 收紧大标题与第一行分图标题之间的距离
            row_gap       = 0.105;
            col_gap       = 0.065;
        else
            left_margin   = 0.090;
            right_margin  = 0.040;
            bottom_margin = 0.090;
            top_margin    = 0.122;
            row_gap       = 0.090;
            col_gap       = 0.060;
        end

        axes_w = (1 - left_margin - right_margin - (ncol-1)*col_gap) / ncol;
        axes_h = (1 - top_margin - bottom_margin - (nrow-1)*row_gap) / nrow;
        if axes_w <= 0 || axes_h <= 0
            error('generate_paper_figures:InvalidLayout', ...
                  'Figure layout margins/gaps leave no room for axes.');
        end

        layout = struct('fig', fig, 'nrow', nrow, 'ncol', ncol, ...
                        'left', left_margin, 'bottom', bottom_margin, ...
                        'top', top_margin, 'row_gap', row_gap, ...
                        'col_gap', col_gap, 'axes_w', axes_w, 'axes_h', axes_h);

        % 整图大标题独立于坐标轴布局，避免与第一行分图标题竞争同一空间。
        if SHOW_MAIN_TITLE && nargin >= 5 && ~isempty(main_title)
            % 大标题单独占据最上方区域，但与下方第一行分图标题仅保留
            % 约 1 个字符的间距，避免显得过空。
            title_gap = 0.008;
            title_box_bottom = 1 - top_margin + title_gap;
            title_box_height = max(0.028, top_margin - title_gap - 0.010);
            annotation(fig, 'textbox', [left_margin, title_box_bottom, ...
                       1-left_margin-right_margin, title_box_height], ...
                       'String', main_title, ...
                       'FontSize', 13, 'FontWeight', 'bold', ...
                       'HorizontalAlignment', 'center', ...
                       'VerticalAlignment', 'middle', ...
                       'EdgeColor', 'none', ...
                       'Interpreter', 'none', ...
                       'FitBoxToText', 'off');
        end
    end

    function ax = next_paper_tile(layout, tile_idx)
        % tile_idx 按“从左到右、从上到下”的顺序编号，与 nexttile 习惯一致。
        row = ceil(tile_idx / layout.ncol);
        col = tile_idx - (row-1)*layout.ncol;

        left = layout.left + (col-1)*(layout.axes_w + layout.col_gap);
        top_edge = 1 - layout.top - (row-1)*(layout.axes_h + layout.row_gap);
        bottom = top_edge - layout.axes_h;

        ax = axes('Parent', layout.fig, 'Units', 'normalized', ...
                  'Position', [left, bottom, layout.axes_w, layout.axes_h]);
        set(ax, 'FontSize', 11, 'Box', 'on', 'Layer', 'top', ...
                'TickDir', 'out', 'XTickLabelMode', 'auto');
        % 固定真正的绘图区位置：即使不同子图的 y 轴刻度字符串长度不同，
        % 三幅图的坐标轴绘图区仍保持严格水平对齐。
        try
            ax.PositionConstraint = 'innerposition';
        catch
            try
                ax.ActivePositionProperty = 'position'; %#ok<SETPROP>
            catch
            end
        end
        grid(ax, 'on');
    end

    function panel_title(ax, txt)
        % 分图标题采用 MATLAB 原生 title，位于分图正上方。
        % 行间距在 create_paper_layout 中显式预留，因此不会与上一幅图的
        % x 轴刻度数字重叠，也不会侵入整图大标题区域。
        % 在所有绘图命令完成后再次强制开启主网格，避免高层绘图函数重置网格。
        grid(ax, 'on');
        set(ax, 'XGrid', 'on', 'YGrid', 'on');
        h = title(ax, txt, ...
                  'FontSize', 11, ...
                  'FontWeight', 'bold', ...
                  'Interpreter', 'none');
        try
            h.Units = 'normalized';
            pos = h.Position;
            pos(2) = 1.018;
            h.Position = pos;
        catch
        end
    end

    function save_data(data_struct, name)
        % MAT 保存全部字段；CSV 仅保存具有统一样本长度的时序数据。
        % 短向量/标量元数据（如 event_centers、event_sigma）单独保存为 metadata CSV，
        % 避免不同变量行数不一致导致 table 构造失败。

        save(fullfile(dat_dir, [name, '_data.mat']), '-struct', 'data_struct');
        fprintf('   ✅ 图 %s 数据保存（.mat）\n', name);

        try
            fields = fieldnames(data_struct);

            % ---------- 1) 自动确定主时序长度 ----------
            sample_lengths = [];
            for i = 1:length(fields)
                val = data_struct.(fields{i});
                if isnumeric(val) && ~isempty(val) && ndims(val) <= 2
                    if isvector(val)
                        sample_lengths(end+1) = numel(val); %#ok<AGROW>
                    else
                        sample_lengths(end+1) = size(val,1); %#ok<AGROW>
                    end
                end
            end

            if isempty(sample_lengths)
                target_rows = 0;
            else
                target_rows = max(sample_lengths);
            end

            % ---------- 2) 时序数据 CSV ----------
            T = table();

            % 元数据采用长表保存：Field / Index / Value
            meta_field = {};
            meta_index = [];
            meta_value = [];

            for i = 1:length(fields)
                field_name = fields{i};
                val = data_struct.(field_name);

                if ~isnumeric(val) || isempty(val) || ndims(val) > 2
                    continue;
                end

                if isvector(val)
                    this_len = numel(val);

                    if this_len == target_rows && target_rows > 1
                        % 无论原来是行向量还是列向量，统一转成列
                        T.(field_name) = val(:);
                    else
                        % 标量/短向量作为元数据保存
                        vals = val(:);
                        for k = 1:numel(vals)
                            meta_field{end+1,1} = field_name; %#ok<AGROW>
                            meta_index(end+1,1) = k; %#ok<AGROW>
                            meta_value(end+1,1) = vals(k); %#ok<AGROW>
                        end
                    end

                else
                    this_len = size(val,1);

                    if this_len == target_rows && target_rows > 0
                        for j = 1:size(val,2)
                            col_name = sprintf('%s_%d', field_name, j);
                            T.(col_name) = val(:,j);
                        end
                    else
                        % 与主时序长度不同的二维数值数据转存为元数据长表
                        vals = val(:);
                        for k = 1:numel(vals)
                            meta_field{end+1,1} = field_name; %#ok<AGROW>
                            meta_index(end+1,1) = k; %#ok<AGROW>
                            meta_value(end+1,1) = vals(k); %#ok<AGROW>
                        end
                    end
                end
            end

            if width(T) > 0
                csv_filename = fullfile(dat_dir, [name, '_data.csv']);
                writetable(T, csv_filename);
                fprintf('   ✅ 图 %s 时序数据保存（.csv）\n', name);
            end

            % ---------- 3) 元数据 CSV ----------
            if ~isempty(meta_field)
                T_meta = table(string(meta_field), meta_index, meta_value, ...
                    'VariableNames', {'Field','Index','Value'});
                meta_filename = fullfile(dat_dir, [name, '_metadata.csv']);
                writetable(T_meta, meta_filename);
                fprintf('   ✅ 图 %s 元数据保存（metadata.csv）\n', name);
            end

        catch ME
            % MATLAB warning 的第一个参数必须是消息ID或字符串，不能直接传 MException。
            warning('generate_paper_figures:CSVExportFailed', ...
                'CSV导出失败（不影响MAT保存）：%s', ME.message);
        end
    end

    % 计算误差
    err_fixed = q_fixed - q_ref;
    err_zhang = q_zhang - q_ref;
    err_cao   = q_cao - q_ref;
    err_adapt = q_adapt - q_ref;
    delta_q_adp  = q_adapt - q_ref;
    delta_dq_adp = dq_adapt - dq_ref;

    Ba_zhang = Ba_zhang(:,1:min(3,size(Ba_zhang,2)));
    Ka_zhang = Ka_zhang(:,1:min(3,size(Ka_zhang,2)));
    Ga_zhang = Ga_zhang(:,1:min(3,size(Ga_zhang,2)));
    Ba_cao = Ba_cao(:,1:min(3,size(Ba_cao,2)));
    Ka_cao = Ka_cao(:,1:min(3,size(Ka_cao,2)));
    Ga_cao = Ga_cao(:,1:min(3,size(Ga_cao,2)));
    Ba_adp = Ba_adp(:,1:min(3,size(Ba_adp,2)));
    Ka_adp = Ka_adp(:,1:min(3,size(Ka_adp,2)));
    Ga_adp = Ga_adp(:,1:min(3,size(Ga_adp,2)));

    %% ---------- 3. Paper Fig.5：重建参考轨迹 ----------
    fprintf('\n📊 Paper Fig.5：重建参考轨迹...\n');

    % Display only one representative rehabilitation cycle.
    cycle_duration_fig5 = sum(params.traj_times);
    idx_fig5_cycle = t_all <= cycle_duration_fig5 + params.dt/2;
    t_fig5 = t_all(idx_fig5_cycle);
    q_ref_fig5 = q_ref(idx_fig5_cycle,:);
    dq_ref_fig5 = dq_ref(idx_fig5_cycle,:);
    ddq_ref_fig5 = ddq_ref(idx_fig5_cycle,:);
    [fig5, tl5] = create_paper_layout('Fig. 5', [100 40 1050 1180], 3, 1, ...
                                     'Reconstructed Reference Trajectories');

    axs5 = gobjects(3,1); lgds5 = gobjects(3,1);
    ax = next_paper_tile(tl5, 1); axs5(1)=ax;
    plot(ax, t_fig5, q_ref_fig5(:,1)*180/pi, '-',  'Color',color_ref,   'LineWidth',1.6, 'DisplayName','Hip'); hold(ax,'on');
    plot(ax, t_fig5, q_ref_fig5(:,2)*180/pi, '--', 'Color',color_fixed, 'LineWidth',1.6, 'DisplayName','Knee');
    plot(ax, t_fig5, q_ref_fig5(:,3)*180/pi, ':',  'Color',color_adapt, 'LineWidth',1.6, 'DisplayName','Ankle');
    ylabel(ax, 'Position (deg)');
    lgds5(1)=legend(ax, 'Location','northeast','AutoUpdate','off');
    panel_title(ax, '(a) Joint Position');

    ax = next_paper_tile(tl5, 2); axs5(2)=ax;
    plot(ax, t_fig5, dq_ref_fig5(:,1)*180/pi, '-',  'Color',color_ref,   'LineWidth',1.6, 'DisplayName','Hip'); hold(ax,'on');
    plot(ax, t_fig5, dq_ref_fig5(:,2)*180/pi, '--', 'Color',color_fixed, 'LineWidth',1.6, 'DisplayName','Knee');
    plot(ax, t_fig5, dq_ref_fig5(:,3)*180/pi, ':',  'Color',color_adapt, 'LineWidth',1.6, 'DisplayName','Ankle');
    ylabel(ax, 'Velocity (deg/s)');
    lgds5(2)=legend(ax, 'Location','northeast','AutoUpdate','off');
    panel_title(ax, '(b) Joint Velocity');

    ax = next_paper_tile(tl5, 3); axs5(3)=ax;
    plot(ax, t_fig5, ddq_ref_fig5(:,1)*180/pi, '-',  'Color',color_ref,   'LineWidth',1.6, 'DisplayName','Hip'); hold(ax,'on');
    plot(ax, t_fig5, ddq_ref_fig5(:,2)*180/pi, '--', 'Color',color_fixed, 'LineWidth',1.6, 'DisplayName','Knee');
    plot(ax, t_fig5, ddq_ref_fig5(:,3)*180/pi, ':',  'Color',color_adapt, 'LineWidth',1.6, 'DisplayName','Ankle');
    ylabel(ax, 'Acceleration (deg/s^2)');
    xlabel(ax, 'Time (s)');
    lgds5(3)=legend(ax, 'Location','northeast','AutoUpdate','off');
    panel_title(ax, '(c) Joint Acceleration');
    paper_figure_utils('align_legends',axs5,lgds5);

    export_fig(fig5,'Fig05_reference_trajectories');
    save_data(struct('t',t_all,'q_ref',q_ref,'dq_ref',dq_ref,'ddq_ref',ddq_ref), ...
              'Fig05_reference_trajectories');

    %% ---------- 4. Paper Fig.6：两种文献比较器与 Proposed 参数演化 ----------
    fprintf('📊 Paper Fig.6：Zhang-2018 / Cao-2024-derived / Proposed 参数演化...\n');
    [fig6, tl6] = create_paper_layout('Fig. 6', [40 25 1550 1250], 3, 3, ...
                                     'Online Admittance-Parameter Evolution');

    ctrl_names = {'Zhang-2018','Cao-2024-derived','Proposed'};
    ctrl_vars = {
        Ba_zhang, Ba_cao, Ba_adp;
        Ka_zhang, Ka_cao, Ka_adp;
        Ga_zhang, Ga_cao, Ga_adp};
    ylabels = {'B_a (N m s/rad)','K_a (N m/rad)','G_a (-)'};
    row_names = {'Virtual Damping B_a','Virtual Stiffness K_a','Interaction-Mapping Gain G_a'};
    joint_colors = [0.0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];
    joint_styles = {'-','--',':'};
    axs6 = gobjects(9,1); lgds6 = gobjects(9,1);
    letters = 'abcdefghi';
    tile_idx = 0;
    cycleT = sum(params.traj_times);

    for rr=1:3
        for cc=1:3
            tile_idx = tile_idx + 1;
            ax=next_paper_tile(tl6,tile_idx); axs6(tile_idx)=ax; hold(ax,'on');
            X = ctrl_vars{rr,cc};
            cycle_axis = t_all / cycleT;
            for j=1:3
                plot(ax,cycle_axis,X(:,j),joint_styles{j},'Color',joint_colors(j,:), ...
                    'LineWidth',1.35,'DisplayName',joints_label{j});
            end

            % G_a can be constant for the literature comparators. MATLAB's
            % automatic y-limits may place the maximum constant trace exactly
            % on the upper axes border, making that line visually disappear.
            % Add a small data-dependent vertical margin to every G_a panel so
            % all Hip/Knee/Ankle traces remain clearly visible.
            if rr==3
                gmin = min(X(:),[],'omitnan');
                gmax = max(X(:),[],'omitnan');
                if isfinite(gmin) && isfinite(gmax)
                    gspan = gmax - gmin;
                    if gspan > 0
                        gpad = max(0.08*gspan, 0.02*max(abs([gmin,gmax])));
                    else
                        gpad = max(0.05*abs(gmax), 0.10);
                    end
                    ylim(ax,[gmin-gpad, gmax+gpad]);
                end
            end

            if cc==2
                for cyc=1:params.loop_num-1
                    xline(ax,cyc,':','HandleVisibility','off','LineWidth',0.75);
                end
            end
            ylabel(ax,ylabels{rr});
            if rr==3, xlabel(ax,'Rehabilitation Cycle'); end
            lgds6(tile_idx)=legend(ax,'Location','northeast','AutoUpdate','off');
            panel_title(ax,sprintf('(%c) %s - %s',letters(tile_idx),ctrl_names{cc},row_names{rr}));
        end
    end
    paper_figure_utils('align_legends',axs6,lgds6);
    export_fig(fig6,'Fig06_admittance_parameters_two_literature_vs_proposed');
    save_data(struct('t',t_all, ...
        'Ba_zhang2018',Ba_zhang,'Ka_zhang2018',Ka_zhang,'Ga_zhang2018',Ga_zhang, ...
        'Ba_cao2024',Ba_cao,'Ka_cao2024',Ka_cao,'Ga_cao2024',Ga_cao, ...
        'Ba_proposed',Ba_adp,'Ka_proposed',Ka_adp,'Ga_proposed',Ga_adp, ...
        'cao_cycle_duration_s',cycleT), ...
        'Fig06_admittance_parameters_two_literature_vs_proposed');

    %% ---------- 5. Paper Fig.7：事件扰动区间轨迹响应 ----------
    fprintf('📊 Paper Fig.7：事件扰动区间轨迹响应...\n');
    [fig7, tl7] = create_paper_layout('Fig. 7', [100 35 1050 1220], 3, 1, ...
                                     'Joint Trajectory Responses during Event Disturbances');

    idx = event_timing.event_mask;
    tw = t_all(idx);
    panel_names = {'(a) Hip','(b) Knee','(c) Ankle'};
    axs7 = gobjects(3,1); lgds7 = gobjects(3,1);
    for j = 1:3
        ax = next_paper_tile(tl7, j); axs7(j)=ax;
        plot(ax, tw, q_ref(idx,j)*180/pi, '-', 'Color',color_ref, 'LineWidth',2.0,'DisplayName','Reference'); hold(ax,'on');
        plot(ax, tw, q_fixed(idx,j)*180/pi, '--','Color',color_fixed,'LineWidth',1.35,'DisplayName','Fixed');
        plot(ax, tw, q_zhang(idx,j)*180/pi, '-.','Color',color_zhang,'LineWidth',1.35,'DisplayName','Zhang-2018 VA');
        plot(ax, tw, q_cao(idx,j)*180/pi, ':','Color',color_cao,'LineWidth',1.55,'DisplayName','Cao-2024-derived VA');
        plot(ax, tw, q_adapt(idx,j)*180/pi, '-', 'Color',color_adapt,'LineWidth',1.35,'DisplayName','Proposed');
        ylabel(ax, [joints_label{j},' Angle (deg)']);
        lgds7(j)=legend(ax,'Location','northeast','AutoUpdate','off');
        panel_title(ax,panel_names{j});
        if j == 3, xlabel(ax,'Time (s)'); end
    end
    paper_figure_utils('align_legends',axs7,lgds7,'NumColumns',2);

    export_fig(fig7,'Fig07_joint_tracking');
    save_data(struct('t_win',tw,'q_ref',q_ref(idx,:),'q_fixed',q_fixed(idx,:), ...
        'q_zhang2018',q_zhang(idx,:),'q_cao2024',q_cao(idx,:),'q_proposed',q_adapt(idx,:)), ...
        'Fig07_joint_tracking');

    %% ---------- 6. Paper Fig.8：真实事件扰动区间跟踪误差 ----------
    fprintf('📊 Paper Fig.8：真实事件扰动区间误差...\n');
    [fig8, tl8] = create_paper_layout('Fig. 8', [100 35 1050 1220], 3, 1, ...
                                     sprintf('Tracking Error during Event Disturbances [%.2f, %.2f] s', ...
                                             event_timing.abnormal_window(1), event_timing.abnormal_window(2)));

    idx2 = event_timing.event_mask;
    tw2 = t_all(idx2);
    panel_names = {'(a) Hip Tracking Error','(b) Knee Tracking Error','(c) Ankle Tracking Error'};
    axs8 = gobjects(3,1); lgds8 = gobjects(3,1);
    for j = 1:3
        ax = next_paper_tile(tl8, j); axs8(j)=ax;
        plot(ax,tw2,err_fixed(idx2,j)*180/pi,'--','Color',color_fixed,'DisplayName','Fixed'); hold(ax,'on');
        plot(ax,tw2,err_zhang(idx2,j)*180/pi,'-.','Color',color_zhang,'DisplayName','Zhang-2018 VA');
        plot(ax,tw2,err_cao(idx2,j)*180/pi,':','Color',color_cao,'LineWidth',1.55,'DisplayName','Cao-2024-derived VA');
        plot(ax,tw2,err_adapt(idx2,j)*180/pi,'-','Color',color_adapt,'DisplayName','Proposed');
        ylabel(ax,[joints_label{j},' Error (deg)']);
        lgds8(j)=legend(ax,'Location','northeast','AutoUpdate','off');
        panel_title(ax,panel_names{j});
        if j == 3, xlabel(ax,'Time (s)'); end
    end
    paper_figure_utils('align_legends',axs8,lgds8,'NumColumns',2);

    export_fig(fig8,'Fig08_tracking_error');
    save_data(struct('t_win',tw2,'event_centers',event_timing.centers,'event_sigma',event_timing.sigma, ...
        'err_fixed',err_fixed(idx2,:),'err_zhang2018',err_zhang(idx2,:), ...
        'err_cao2024',err_cao(idx2,:),'err_proposed',err_adapt(idx2,:)), ...
        'Fig08_tracking_error');

    %% ---------- 7. Paper Fig.9：运动平滑度（瞬时 jerk） ----------
    fprintf('📊 Paper Fig.9：事件扰动区间 jerk 响应...\n');
    [fig9, tl9] = create_paper_layout('Fig. 9', [100 35 1050 1220], 3, 1, ...
                                     'Instantaneous Jerk during Event Disturbances');

    jerk_fixed_all = zeros(size(ddq_fixed));
    jerk_zhang_all = zeros(size(ddq_zhang));
    jerk_cao_all   = zeros(size(ddq_cao));
    jerk_adapt_all = zeros(size(ddq_adapt));
    for j = 1:3
        jerk_fixed_all(:,j) = gradient(ddq_fixed(:,j),dt);
        jerk_zhang_all(:,j) = gradient(ddq_zhang(:,j),dt);
        jerk_cao_all(:,j)   = gradient(ddq_cao(:,j),dt);
        jerk_adapt_all(:,j) = gradient(ddq_adapt(:,j),dt);
    end

    idx3 = event_timing.event_mask;
    tw3 = t_all(idx3);
    panel_names = {'(a) Hip Jerk','(b) Knee Jerk','(c) Ankle Jerk'};
    axs9 = gobjects(3,1); lgds9 = gobjects(3,1);
    for j = 1:3
        ax = next_paper_tile(tl9,j); axs9(j)=ax;
        plot(ax,tw3,jerk_fixed_all(idx3,j)*180/pi,'--','Color',color_fixed,'LineWidth',1.30,'DisplayName','Fixed'); hold(ax,'on');
        plot(ax,tw3,jerk_zhang_all(idx3,j)*180/pi,'-.','Color',color_zhang,'LineWidth',1.30,'DisplayName','Zhang-2018 VA');
        plot(ax,tw3,jerk_cao_all(idx3,j)*180/pi,':','Color',color_cao,'LineWidth',1.50,'DisplayName','Cao-2024-derived VA');
        plot(ax,tw3,jerk_adapt_all(idx3,j)*180/pi,'-','Color',color_adapt,'LineWidth',1.30,'DisplayName','Proposed');
        ylabel(ax,[joints_label{j},' Jerk (deg/s^3)']);
        lgds9(j)=legend(ax,'Location','northeast','AutoUpdate','off');
        panel_title(ax,panel_names{j});
        if j == 3, xlabel(ax,'Time (s)'); end
    end
    paper_figure_utils('align_legends',axs9,lgds9,'NumColumns',2);

    export_fig(fig9,'Fig09_motion_smoothness');
    save_data(struct('t_win',tw3,'event_centers',event_timing.centers,'event_sigma',event_timing.sigma, ...
        'jerk_fixed',jerk_fixed_all(idx3,:),'jerk_zhang2018',jerk_zhang_all(idx3,:), ...
        'jerk_cao2024',jerk_cao_all(idx3,:),'jerk_proposed',jerk_adapt_all(idx3,:), ...
        'note','Jerk was computed on the full time series before event-window extraction.'), ...
        'Fig09_motion_smoothness');

    %% ---------- 8. Paper Fig.10：三项 primary performance metrics ----------
    fprintf('📊 Paper Fig.10：四控制器三项 primary performance metrics...\n');
    [fig10, tl10] = create_paper_layout('Fig. 10', [70 90 1450 520], 1, 3, ...
                                       'Primary Quantitative Performance');

    rmse_f = metrics_fixed.RMSE(:); rmse_z = metrics_zhang.RMSE(:); rmse_c = metrics_cao.RMSE(:); rmse_a = metrics_adaptive.RMSE(:);
    jerk_f = metrics_fixed.RMS_jerk_deg(:); jerk_z = metrics_zhang.RMS_jerk_deg(:); jerk_c = metrics_cao.RMS_jerk_deg(:); jerk_a = metrics_adaptive.RMS_jerk_deg(:);
    raw_f = metrics_fixed.RMS_InteractionTorque(:); raw_z = metrics_zhang.RMS_InteractionTorque(:); raw_c = metrics_cao.RMS_InteractionTorque(:); raw_a = metrics_adaptive.RMS_InteractionTorque(:);
    legends4 = {'Fixed','Zhang-2018 VA','Cao-2024-derived VA','Proposed'};

    axs10=gobjects(3,1); lgds10=gobjects(3,1);
    ax = next_paper_tile(tl10,1); axs10(1)=ax;
    bar(ax,[rmse_f,rmse_z,rmse_c,rmse_a],'grouped'); set(ax,'XTick',1:3,'XTickLabel',joints_label);
    ylabel(ax,'RMSE (deg)'); lgds10(1)=legend(ax,legends4,'Location','northeast','AutoUpdate','off'); panel_title(ax,'(a) Tracking RMSE');

    ax = next_paper_tile(tl10,2); axs10(2)=ax;
    bar(ax,[jerk_f,jerk_z,jerk_c,jerk_a],'grouped'); set(ax,'XTick',1:3,'XTickLabel',joints_label);
    ylabel(ax,'Jerk RMS (deg/s^3)'); lgds10(2)=legend(ax,legends4,'Location','northeast','AutoUpdate','off'); panel_title(ax,'(b) Motion Smoothness');

    ax = next_paper_tile(tl10,3); axs10(3)=ax;
    bar(ax,[raw_f,raw_z,raw_c,raw_a],'grouped'); set(ax,'XTick',1:3,'XTickLabel',joints_label);
    ylabel(ax,'Raw Residual RMS (N m)'); lgds10(3)=legend(ax,legends4,'Location','northeast','AutoUpdate','off'); panel_title(ax,'(c) Interaction Residual');
    paper_figure_utils('align_legends',axs10,lgds10);

    export_fig(fig10,'Fig10_quantitative_performance_four_controller');
    fig10data = struct('rmse_fixed',rmse_f,'rmse_zhang2018',rmse_z,'rmse_cao2024',rmse_c,'rmse_proposed',rmse_a, ...
        'jerk_fixed_deg_s3',jerk_f,'jerk_zhang2018_deg_s3',jerk_z,'jerk_cao2024_deg_s3',jerk_c,'jerk_proposed_deg_s3',jerk_a, ...
        'raw_fixed',raw_f,'raw_zhang2018',raw_z,'raw_cao2024',raw_c,'raw_proposed',raw_a, ...
        'note','Primary performance contains RMSE, jerk RMS, and raw interaction-residual RMS only. Position saturation is reported separately in Fig. S1.');
    save_data(fig10data,'Fig10_quantitative_performance_four_controller');

    %% ---------- 9. Paper Fig.11：四控制器 Monte Carlo 鲁棒性 ----------
    fprintf('📊 Paper Fig.11：四控制器 Monte Carlo 鲁棒性箱线图...\n');
    if ~isempty(monte_carlo_path) && exist(monte_carlo_path,'dir')
        mc_mat = dir(fullfile(monte_carlo_path,'MC_stats_*.mat'));
        if ~isempty(mc_mat)
            mc_data = load(fullfile(mc_mat(1).folder,mc_mat(1).name));
            if isfield(mc_data,'stats') && isfield(mc_data.stats,'fixed') && ...
                    isfield(mc_data.stats,'zhang2018') && isfield(mc_data.stats,'cao2024') && ...
                    isfield(mc_data.stats,'adaptive')
                st = mc_data.stats;
                if isfield(st,'success_mask')
                    valid_mc = logical(st.success_mask(:));
                else
                    valid_mc = all(isfinite(st.fixed.RMSE),2) & ...
                               all(isfinite(st.zhang2018.RMSE),2) & ...
                               all(isfinite(st.cao2024.RMSE),2) & ...
                               all(isfinite(st.adaptive.RMSE),2);
                end
                n_valid_mc = sum(valid_mc);
                if n_valid_mc < 2
                    warning('MC stats 仅有 %d 个有效四控制器配对 trial，跳过 Fig.11',n_valid_mc);
                else
                    [fig11, tl11] = create_paper_layout('Fig. 11',[45 20 1550 1250],3,3, ...
                        sprintf('Paired Monte Carlo Robustness (n=%d valid trials)',n_valid_mc));
                    letters='abcdefghi'; panel_idx=0; labels={'Fixed','Zhang18','Cao24','Proposed'};
                    for j=1:3
                        panel_idx=panel_idx+1; ax=next_paper_tile(tl11,panel_idx);
                        f=st.fixed.RMSE(valid_mc,j); z=st.zhang2018.RMSE(valid_mc,j); c=st.cao2024.RMSE(valid_mc,j); a=st.adaptive.RMSE(valid_mc,j);
                        boxplot(ax,[f;z;c;a],[zeros(size(f));ones(size(z));2*ones(size(c));3*ones(size(a))],'Labels',labels);
                        ylabel(ax,'RMSE (deg)'); panel_title(ax,sprintf('(%c) %s RMSE',letters(panel_idx),joints_label{j}));
                    end
                    for j=1:3
                        panel_idx=panel_idx+1; ax=next_paper_tile(tl11,panel_idx);
                        f=st.fixed.JerkRMS(valid_mc,j)*180/pi; z=st.zhang2018.JerkRMS(valid_mc,j)*180/pi; c=st.cao2024.JerkRMS(valid_mc,j)*180/pi; a=st.adaptive.JerkRMS(valid_mc,j)*180/pi;
                        boxplot(ax,[f;z;c;a],[zeros(size(f));ones(size(z));2*ones(size(c));3*ones(size(a))],'Labels',labels);
                        ylabel(ax,'Jerk RMS (deg/s^3)'); panel_title(ax,sprintf('(%c) %s Jerk RMS',letters(panel_idx),joints_label{j}));
                    end
                    for j=1:3
                        panel_idx=panel_idx+1; ax=next_paper_tile(tl11,panel_idx);
                        f=st.fixed.RMS_InteractionTorque(valid_mc,j); z=st.zhang2018.RMS_InteractionTorque(valid_mc,j); c=st.cao2024.RMS_InteractionTorque(valid_mc,j); a=st.adaptive.RMS_InteractionTorque(valid_mc,j);
                        boxplot(ax,[f;z;c;a],[zeros(size(f));ones(size(z));2*ones(size(c));3*ones(size(a))],'Labels',labels);
                        ylabel(ax,'Raw Residual RMS (N m)'); panel_title(ax,sprintf('(%c) %s Residual',letters(panel_idx),joints_label{j}));
                    end
                    export_fig(fig11,'Fig11_monte_carlo_robustness_four_controller');

                    fig11_data=struct('valid_trial_indices',find(valid_mc),'n_valid_trials',n_valid_mc);
                    for j=1:3
                        nm=joints_label{j};
                        fig11_data.(['fixed_rmse_' nm])=st.fixed.RMSE(valid_mc,j);
                        fig11_data.(['zhang2018_rmse_' nm])=st.zhang2018.RMSE(valid_mc,j);
                        fig11_data.(['cao2024_rmse_' nm])=st.cao2024.RMSE(valid_mc,j);
                        fig11_data.(['proposed_rmse_' nm])=st.adaptive.RMSE(valid_mc,j);
                        fig11_data.(['fixed_jerk_' nm])=st.fixed.JerkRMS(valid_mc,j)*180/pi;
                        fig11_data.(['zhang2018_jerk_' nm])=st.zhang2018.JerkRMS(valid_mc,j)*180/pi;
                        fig11_data.(['cao2024_jerk_' nm])=st.cao2024.JerkRMS(valid_mc,j)*180/pi;
                        fig11_data.(['proposed_jerk_' nm])=st.adaptive.JerkRMS(valid_mc,j)*180/pi;
                        fig11_data.(['fixed_raw_' nm])=st.fixed.RMS_InteractionTorque(valid_mc,j);
                        fig11_data.(['zhang2018_raw_' nm])=st.zhang2018.RMS_InteractionTorque(valid_mc,j);
                        fig11_data.(['cao2024_raw_' nm])=st.cao2024.RMS_InteractionTorque(valid_mc,j);
                        fig11_data.(['proposed_raw_' nm])=st.adaptive.RMS_InteractionTorque(valid_mc,j);
                    end
                    save_data(fig11_data,'Fig11_monte_carlo_robustness_four_controller');

                    if isfield(st.fixed,'PositionSaturationPct') && ...
                            isfield(st.zhang2018,'PositionSaturationPct') && ...
                            isfield(st.cao2024,'PositionSaturationPct') && ...
                            isfield(st.adaptive,'PositionSaturationPct')
                        [figS1, tlS1] = create_paper_layout('Fig. S1',[70 80 1450 500],1,3, ...
                            sprintf('Supervisory Diagnostic: Position-Correction Saturation (n=%d)',n_valid_mc));
                        sat_data=struct('valid_trial_indices',find(valid_mc),'n_valid_trials',n_valid_mc);
                        for j=1:3
                            ax=next_paper_tile(tlS1,j);
                            f=st.fixed.PositionSaturationPct(valid_mc,j);
                            z=st.zhang2018.PositionSaturationPct(valid_mc,j);
                            c=st.cao2024.PositionSaturationPct(valid_mc,j);
                            a=st.adaptive.PositionSaturationPct(valid_mc,j);
                            boxplot(ax,[f;z;c;a],[zeros(size(f));ones(size(z));2*ones(size(c));3*ones(size(a))],'Labels',labels);
                            ylabel(ax,'Position saturation (%)'); panel_title(ax,sprintf('(%c) %s',char('a'+j-1),joints_label{j}));
                            nm=joints_label{j};
                            sat_data.(['fixed_sat_' nm])=f;
                            sat_data.(['zhang2018_sat_' nm])=z;
                            sat_data.(['cao2024_sat_' nm])=c;
                            sat_data.(['proposed_sat_' nm])=a;
                        end
                        export_fig(figS1,'FigS1_position_saturation_diagnostic');
                        save_data(sat_data,'FigS1_position_saturation_diagnostic');
                    end
                end
            else
                warning('MC stats 缺少 fixed/zhang2018/cao2024/adaptive 字段，跳过 Fig.11');
            end
        else
            warning('MC_stats 文件未找到，跳过 Fig.11');
        end
    else
        fprintf('未找到蒙特卡洛结果，跳过 Fig.11\n');
    end

    fprintf('\n🎉 论文图表全部生成完成！\n');
end