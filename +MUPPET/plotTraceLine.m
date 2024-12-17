function ax = plotTraceLine(varargin)
% plots the trace line of a call against a spectrogram.
%
% Written by Wilfried Beslin
% Last updated 2024-12-17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check if first argument is an Axes object or not
    arg1 = varargin{1};
    if isa(arg1, 'matlab.graphics.axis.Axes')
        ax = arg1;
        plot_args = varargin(2:end);
    else
        ax = [];
        plot_args = varargin;
    end

    % parse remaining arguments
    p = inputParser();
    p.addRequired('t_stft', @(v)validateattributes(v,{'numeric'},{'row'}))
    p.addRequired('f_stft', @(v)validateattributes(v,{'numeric'},{'column'}))
    p.addRequired('psdm', @(v)validateattributes(v,{'numeric'},{'2d'}))
    p.addRequired('t_trace', @(v)validateattributes(v,{'cell','numeric'},{})) % may be a row vector, or a cell of row vecrtors where each cell element corresponds to a seperate trace line.
    p.addRequired('f_trace', @(v)validateattributes(v,{'cell','numeric'},{})) % may be a column vector, or a cell of column vecrtors where each cell element corresponds to a seperate trace line.
    p.addOptional('annot_f_range', [], @(v)validateattributes(v,{'numeric'},{'numel',2,'increasing'}))
    p.addParameter('LineData', struct.empty(0,0), @isstruct) % struct array containing line properties, e.g. LineWidth (array has 1 element per trace)
    p.addParameter('CAxis', [], @(v)validateattributes(v,{'numeric'},{'increasing','numel',2})) % determines the bounds of the colour scale
    p.addParameter('LogCols', true, @(v)validateattributes(v,{'logical'},{'scalar'})) % determines if colour scale is logged or not
    
    p.parse(plot_args{:});
    t_stft = p.Results.t_stft;
    f_stft = p.Results.f_stft;
    psdm = p.Results.psdm;
    t_trace = p.Results.t_trace;
    if isnumeric(t_trace)
        t_trace = {t_trace};
    end
    f_trace = p.Results.f_trace;
    if isnumeric(f_trace)
        f_trace = {f_trace};
    end
    
    num_traces = numel(f_trace);
    annot_f_range = p.Results.annot_f_range;
    line_data = p.Results.LineData;
    caxis_val = p.Results.CAxis;
    do_logcols = p.Results.LogCols;
    
    % get a suitable maximum value to plot elements above the spectrogram
    z_top = ceil(max(psdm(:)));
    
    % adjust data variables to generate an aesthetic spectrogram plot
    dt = mean(diff(t_stft));
    t_surf = t_stft - dt/2;
    t_surf = [t_surf, t_surf(end) + dt];
    
    df = mean(diff(f_stft));
    f_surf = f_stft - df/2;
    f_surf = [f_surf; f_surf(end) + df];
    
    psdm_surf = psdm;
    psdm_surf = [psdm_surf, psdm_surf(:,end)];
    psdm_surf = [psdm_surf; psdm_surf(end,:)];
    
    % initialize axes
    if isempty(ax)
        ax = axes();
    end
    ax.NextPlot = 'add';
    
    % plot the spectrogram
    surf(ax, t_surf, f_surf, psdm_surf, 'EdgeColor','none');
    axis(ax, 'xy');
    view(ax, 0, 90);
    if ~isempty(caxis_val)
        caxis(ax, caxis_val)
    end
    if do_logcols
        set(ax, 'ColorScale', 'log')
    end
    
    % add frequency bounds if available
    if ~isempty(annot_f_range)
        plot3(ax, t_stft([1,end]), repelem(annot_f_range(1),1,2), [z_top,z_top], 'w:', 'LineWidth',1.5)
        plot3(ax, t_stft([1,end]), repelem(annot_f_range(2),1,2), [z_top,z_top], 'w:', 'LineWidth',1.5)
    end
    
    % plot trace line(s)
    line_data_fields = fieldnames(line_data);
    num_line_data_fields = numel(line_data_fields);
    line_obj = cell(1,num_traces);
    for ii = 1:num_traces
        t_trace_ii = t_trace{ii};
        f_trace_ii = f_trace{ii};
        
        if isempty(line_data)
            line_vars = {'o-', 'LineWidth',1, 'MarkerSize',4};
        else
            line_vars = cell(1, num_line_data_fields*2);
            line_vars(1:2:end) = line_data_fields;
            line_vars(2:2:end) = struct2cell(line_data(ii));
        end
       line_obj{ii} =  plot3(ax, t_trace_ii, f_trace_ii', z_top.*ones(size(t_trace_ii)), line_vars{:});
    end
    line_obj = [line_obj{:}];
    
    % create legend if there is more than one trace
    if num_traces > 1
        lgd = legend(ax, line_obj, 'Location','northeast');
        lgd.AutoUpdate = 'off';
    end
    
    %** TEST
    % create lowess-smoothed versions of trace and plot them
    %f_trace_smooth_narrow = smoothdata(f_trace, 'lowess', 7);
    %f_trace_smooth_wide = smoothdata(f_trace, 'loess', 15);
    %plot3(ax, t_trace, f_trace_smooth_narrow', ones(size(t_trace)), 'kx-', 'LineWidth',1, 'MarkerSize',4)
    %plot3(ax, t_trace, f_trace_smooth_wide', ones(size(t_trace)), 'wx-', 'LineWidth',1, 'MarkerSize',4)
    %** END TEST
    
    % set axis labels
    xlabel(ax, 'Time [s]')
    ylabel(ax, 'Frequency [Hz]')
end