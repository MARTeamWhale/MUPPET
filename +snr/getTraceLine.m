function i_f_trace = getTraceLine(f_stft, psdm, penalty_coeff)
% Find an ideal trace line through a call in a spectrogram by using a
% "shortest path" searching algorithm (Dijkstra's algorithm)
%
% Written by Wilfried Beslin
% Last updated 2024-07-15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - consider adding trace smoothing options as optional input arguments.
% But know that in this case, can no longer return frequency indices;
% rather, would need to return both interpolated frequencies and
% magnitudes.

    % define the penalty function for jumps across fequency
    penalty_fcn = @(df) penalty_coeff.*(df.^2) + 1;
    
    % get counts from spectrogram
    [nf, nt] = size(psdm);

    % find the peak of the spectrogram
    [~, i_peak] = max(psdm(:));
    [r_peak, c_peak] = ind2sub(size(psdm), i_peak);
    
    % create a matrix of initial node weights based on PSD intensity
    all_og_weights = 10*log10(psdm);
    all_og_weights = abs(all_og_weights - all_og_weights(i_peak));
    
    % create directional networks representing all possible paths from the
    % peak to the end (forward) and beginning (backward) of the spectrogram
    
    %%% forward
    fw_og_weights = all_og_weights(:,(c_peak+1):end);
    fw_net = makePathNetwork(r_peak, fw_og_weights, f_stft, penalty_fcn);
    
    %%% backward
    bw_og_weights = all_og_weights(:,(c_peak-1):-1:1);
    bw_net = makePathNetwork(r_peak, bw_og_weights, f_stft, penalty_fcn);
    
    % find the path in the forward graph with smallest total weight
    fw_trace_nodes = getShortestPath(fw_net, nf);
    
    % do the same for the backwards network
    bw_trace_nodes = getShortestPath(bw_net, nf);
    
    % translate the trace node IDs into time-frequency points
    % (only the frequency indices are returned)
    [r_trace_fw, c_trace_fw] = ind2sub([nf,nt-c_peak], fw_trace_nodes(2:end) - 1);
    [r_trace_bw, c_trace_bw] = ind2sub([nf,c_peak-1], bw_trace_nodes(2:end) - 1);
    
    i_f_trace = [fliplr(r_trace_bw), r_peak, r_trace_fw];
    %i_t_trace = [c_trace_bw, c_peak, c_trace_fw+c_peak];
    
    %t_trace = t_stft(i_t_trace);
    %f_trace = f_stft(i_f_trace);
end


%% makePathNetwork --------------------------------------------------------
function path_net = makePathNetwork(r_start, jumpnode_weights, penalty_var, penalty_fcn)
% Create a weighted directional graph representing every possible path from
% a starting point through to the end of a rectangular network of nodes.
% Intended to represent (partial) trace lines through a spectrogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Here, "jump node" means any node except the starting node.
    % Jump nodes are arranged in a rectangular matrix, and paths flow
    % through every node from left to right.

    % get counts from the jump node matrix
    [n_rows, n_jumpnode_cols] = size(jumpnode_weights);
    n_jumpnodes = prod([n_rows,n_jumpnode_cols]);
    n_jumpnodes_connected = n_jumpnodes - n_rows;
    
    % create matrix of indices for all jump nodes
    jumpnode_idxmat = reshape(1:n_jumpnodes, n_rows, n_jumpnode_cols) + 1;
    
    % create matrices representing jumps from source nodes to target nodes.
    % Rows = target nodes, Columns = source nodes.
    
    %%% source
    jumpnode_src_mat = repmat((1:n_jumpnodes_connected)+1, n_rows, 1);
    full_src_mat = [ones(n_rows,1), jumpnode_src_mat];
    
    %%% target
    jumpnode_tgt_mat = repmat((0:(n_rows-1))', 1, n_jumpnodes_connected) + n_rows.*repelem(1:(n_jumpnode_cols-1), 1, n_rows) + 2;
    full_tgt_mat = [jumpnode_idxmat(:,1), jumpnode_tgt_mat];
    
    % create matrix of original weights for each jump
    % (this consists simply of the weight of the corresponding target node)
    jumpnode_ogweight_mat = repelem(jumpnode_weights(:,2:end), 1, n_rows);
    full_ogweight_mat = [jumpnode_weights(:,1), jumpnode_ogweight_mat];
    
    % create matrix of penalty factors for each jump
    ogpenalty_refmat = abs(repmat(penalty_var, 1, n_rows) - repmat(penalty_var', n_rows, 1));
    jumpnode_ogpenalty_mat = repmat(ogpenalty_refmat, 1, n_jumpnode_cols-1);
    full_ogpenalty_mat = [abs(penalty_var - penalty_var(r_start)), jumpnode_ogpenalty_mat];
    full_penalty_mat = penalty_fcn(full_ogpenalty_mat);
    
    % combine the original weight matrix with the penalty matrix to create
    % the final edge weights
    full_weight_mat = full_ogweight_mat .* full_penalty_mat;
    
    % create the source node, target node, and weight vectors that will be
    % used to construct the directional graph object
    src_vec = reshape(full_src_mat, 1, numel(full_src_mat));
    tgt_vec = reshape(full_tgt_mat, 1, numel(full_tgt_mat));
    weight_vec = reshape(full_weight_mat, 1, numel(full_weight_mat));
    
    % create digraph object
    path_net = digraph(src_vec, tgt_vec, weight_vec);
end


%% getShortestPath --------------------------------------------------------
function shortest_path = getShortestPath(path_net, nf)
% Return the "shortest" path from a source node to the edge of a
% rectangular flow network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_nodes = height(path_net.Nodes);

    % get target nodes
    tgt = ((n_nodes-nf)+1):n_nodes;
    
    % get shortest paths to all target nodes
    [paths, dists] = shortestpathtree(path_net, 1, tgt, 'OutputForm','cell');
    
    % identify the shortest of all paths
    [~, i_shortest] = min(dists);
    shortest_path = paths{i_shortest};
end