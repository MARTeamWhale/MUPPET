function ax = plotTraceLine(varargin)
% plots the trace line of a call against a spectrogram.
%
% Written by Wilfried Beslin
% Last updated 2024-07-26
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
    p.addRequired('t_trace', @(v)validateattributes(v,{'numeric'},{'row'}))
    p.addRequired('f_trace', @(v)validateattributes(v,{'numeric'},{'column'}))
    p.addOptional('annot_f_range', [], @(v)validateattributes(v,{'numeric'},{'numel',2,'increasing'}))
    
    p.parse(plot_args{:});
    t_stft = p.Results.t_stft;
    f_stft = p.Results.f_stft;
    psdm = p.Results.psdm;
    t_trace = p.Results.t_trace;
    f_trace = p.Results.f_trace;
    annot_f_range = p.Results.annot_f_range;
    
    % transform the spectrogram to be in log scale, with a max of 0 dB
    psdm_log = 10.*log10(psdm);
    psdm_log = psdm_log - max(psdm_log(:));
    
    % adjust data variables to generate an aesthetic spectrogram plot
    dt = mean(diff(t_stft));
    t_surf = t_stft - dt/2;
    t_surf = [t_surf, t_surf(end) + dt];
    
    df = mean(diff(f_stft));
    f_surf = f_stft - df/2;
    f_surf = [f_surf; f_surf(end) + df];
    
    psdm_surf = psdm_log;
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
    caxis(ax, [-60, -3])
    set(ax, 'ColorScale', 'log')
    
    % add frequency bounds if available
    if ~isempty(annot_f_range)
        plot3(ax, t_stft([1,end]), repelem(annot_f_range(1),1,2), [1,1], 'w:', 'LineWidth',1.5)
        plot3(ax, t_stft([1,end]), repelem(annot_f_range(2),1,2), [1,1], 'w:', 'LineWidth',1.5)
    end
    
    % plot trace line
    plot3(ax, t_trace, f_trace', ones(size(t_trace)), 'ro-', 'LineWidth',1, 'MarkerSize',4)
    
    % set axis labels
    xlabel(ax, 'Time [s]')
    ylabel(ax, 'Frequency [Hz]')
end