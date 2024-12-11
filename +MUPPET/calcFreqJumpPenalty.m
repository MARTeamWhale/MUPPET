function penalty_scores = calcFreqJumpPenalty(df, df_mu, ddf_sigma, penalty_at_sigma, varargin)
% Transforms an input array of frequency differences into penalty scores
% (intended to be in dB), to be used for tracing time-frequency changes
% within a call. 
% 
% This function works by applying a Normal probability distribution to the
% set of frequency differences centred at some value "df_mu", which
% corresponds to the most likely frequency difference (often 0 Hz). The
% probability densities are then converted into decibel penalty scores
% ranging from 0 dB when "df" = "df_mu" and progressing to infinite dB as
% "df" moves away from "df_mu" in either direction. This conversion is done
% by normalizing the probability densities such that they range from 0 to
% 1, and then applying a negative logarithmic function with a coefficient
% determined by the parameter "penalty_at_sigma".
%
% "df" can be a vector or matrix. If a matrix, scores will be calculated
% for each column.
%
% Scores can optionally be plotted, but only if the input is a vector.
%
% Written by Wilfried Beslin
% Last updated 2024-12-11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - If I want to support having a different df_mu for each column in df,
% then I will need to allow df_mu to be a vector or scalar. Will also need
% to revise plotting if df_mu can be a vector.

    % parse df
    row_input = isrow(df);
    if row_input
        df = df';
    end
    
    % parse optional input
    p = inputParser();
    p.addParameter('PlotScores', false, @(v)validateattributes(v,{'logical'},{'scalar'}));
    p.parse(varargin{:})
    do_plot = p.Results.PlotScores;

    % define the Normal distribution (Gaussian) PDF function
    gaussFcn = @(x) (1./(ddf_sigma.*sqrt(2.*pi))).*exp((-(x-df_mu).^2)./(2.*ddf_sigma.^2));

    % get probabilistic scores for the frequency differences using the
    % Normal distribution
    pd_scores = gaussFcn(df);

    % normalize the scores such that they are 1 at df = 0
    pd_scores_norm = pd_scores./max(pd_scores, [], 1);

    % calculate the coefficient for the conversion function
    norm_pd_at_sigma = gaussFcn(df_mu + ddf_sigma)/gaussFcn(df_mu);
    converter_coeff = -penalty_at_sigma./log(norm_pd_at_sigma);

    % convert the probabilistic scores
    penalty_scores = -converter_coeff.*log(pd_scores_norm);

    % plot scores if requested (only for vector input)
    if do_plot
        if isvector(df)
            x_label = '\Deltaf [Hz]';
            y_label = 'Penalty [dB]';

            ax = gca();
            if strcmp(ax.XLabel.String, x_label) && strcmp(ax.YLabel.String, y_label)
                ax.NextPlot = 'add';
            else
                cla(ax);
            end

            plot(df, penalty_scores, 'o-');

            grid(ax, 'on')
            box(ax, 'on')

            xlabel(ax, x_label)
            ylabel(ax, y_label)
        else
            warning('Scores can only be plotted if ''df'' is a vector')
        end
    end

    % transpose output if needed
    if row_input
        penalty_scores = penalty_scores';
    end
end