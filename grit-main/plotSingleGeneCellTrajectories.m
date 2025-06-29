function plotSingleGeneCellTrajectories(scdata, Tgrid, gene_index)
% plotSingleGeneCellTrajectories Visualizes the expression trajectory for a
%   specific gene across a subset of individual cells over time.
%
%   plotSingleGeneCellTrajectories(scdata, Tgrid, gene_index) plots the
%   expression trajectories of the specified 'gene_index' for the first
%   5 cells. If the total number of cells is less than 5, all available
%   cells are plotted. It assumes the number of cells remains constant
%   across all time points.
%
%   Inputs:
%     scdata     - A cell array where each element 'scdata{t}' is a matrix
%                  representing gene expression at time point 't'.
%                  Dimensions: (number_of_genes x number_of_cells).
%                  Assumes 'number_of_cells' is constant across all 't'.
%     Tgrid      - A numeric vector containing the time points.
%                  Length must match the number of elements in 'scdata'.
%     gene_index - An integer specifying the index of the gene to visualize.
%                  Must be between 1 and the total number of genes.
%
%   Example Usage (assuming 'scdata' and 'Tgrid' are from a simulation):
%     % [scdata_example, Tgrid_example] = linData_custom(simOpts, A, b, x0);
%     % plotSingleGeneCellTrajectories(scdata_example, Tgrid_example, 3); % Plots gene 3

    % --- Input Validation ---
    if isempty(scdata) || isempty(Tgrid)
        error('plotSingleGeneCellTrajectories:EmptyInput', ...
              'Input ''scdata'' and ''Tgrid'' cannot be empty.');
    end

    if length(scdata) ~= length(Tgrid)
        error('plotSingleGeneCellTrajectories:MismatchedTimePoints', ...
              'The number of time points in ''scdata'' must match the length of ''Tgrid''.');
    end

    % Determine total number of genes and cells from the first time point
    n_genes = size(scdata{1}, 1);
    n_cells = size(scdata{1}, 2); % Assume number of cells remains constant

    if gene_index < 1 || gene_index > n_genes || mod(gene_index, 1) ~= 0
        error('plotSingleGeneCellTrajectories:InvalidGeneIndex', ...
              '''gene_index'' must be a positive integer between 1 and %d (number of genes).', n_genes);
    end

    n_timepoints = length(scdata);

    % --- Extract Data for the Specified Gene Across All Cells ---
    % Initialize a matrix to store the expression of the chosen gene for all cells
    % across all time points.
    % Rows: Time Points, Columns: Cells
    gene_expression_all_cells = zeros(n_timepoints, n_cells);

    for t = 1:n_timepoints
        current_time_data = scdata{t};
        % Check for empty time point data (shouldn't happen with constant cell assumption,
        % but good for robustness if data could be sparse)
        if ~isempty(current_time_data)
            % Extract the row corresponding to the specified gene_index
            gene_expression_all_cells(t, :) = current_time_data(gene_index, :);
        end
    end

    % --- Determine Number of Cells to Graph ---
    cells_to_graph = min(5, n_cells);

    % --- Plotting Individual Cell Trajectories for the Specified Gene ---
    figure; % Create a new figure window
    hold on; % Allow multiple plots on the same axes
    grid on; % Add a grid to the plot

    % Get a colormap for distinct lines for each cell
    colors = lines(cells_to_graph);

    for cell_idx = 1:cells_to_graph % Iterate through the first 5 cells or all if less than 5
        cell_trajectory = gene_expression_all_cells(:, cell_idx);

        plot(Tgrid, cell_trajectory, ...
             'Color', colors(cell_idx, :), ...
             'LineWidth', 1, ...
             'DisplayName', ['Cell ' num2str(cell_idx)]);
    end

    % --- Customize Plot Appearance ---
    title(['Expression Trajectories for Gene ' num2str(gene_index) ' Across Cells']);
    xlabel('Time');
    ylabel(['Expression of Gene ' num2str(gene_index)]);
    legend('show', 'Location', 'bestoutside'); % Display legend outside plot area
    set(gca, 'FontSize', 12); % Set general font size for axes labels and ticks
    box on; % Display the box around the plot area
    hold off; % Release the plot for further modifications

end
