%************************************************************************
% Gray–Scott model 2D simulation in GNU Octave
% Created by HA2ERZ & ChatGPT in 2025
% Released under CC0 1.0 Universal (Public Domain Dedication)
% You may copy, modify, distribute and perform the work,
% even for commercial purposes, without asking permission.
%************************************************************************

clear; clc;
pkg load image;

% -----------------------------------------------------------------------
% Mathematical background (short version)
% -----------------------------------------------------------------------
% The Gray-Scott reaction-diffusion system models two concentrations, U(x,y,t)
% and V(x,y,t):
%
%   dU/dt = Du * ΔU - U*V^2 + F*(1-U)
%   dV/dt = Dv * ΔV + U*V^2 - (F+k)*V
%
% where:
%   - Du, Dv: diffusion coefficients
%   - F: feed rate of U
%   - k: kill rate of V
%   - Δ: 2D Laplacian operator (here implemented by finite differences)
%
% Time integration below uses explicit Euler with internal substeps
% (dt_sub <= max_dt) for better stability:
%   U_{n+1} = U_n + dt_sub * RHS_U(U_n, V_n)
%   V_{n+1} = V_n + dt_sub * RHS_V(U_n, V_n)
%
% Error propagation is estimated by linearized sensitivity equations:
%   d/dp [U_{n+1};V_{n+1}] ≈ J_state * d/dp [U_n;V_n] + dRHS/dp
% for parameters p in {F, k}. This yields uncertainty maps from user-supplied
% parameter uncertainties sigma_F and sigma_k.

% --- Basic settings ---
n = 300;                % Grid size (tested up to 1200)
Du = 0.16;              % Diffusion coefficient for U
Dv = 0.08;              % Diffusion coefficient for V
dt = 1.0;               % Time step
steps = 5000;           % Number of iterations
save_interval = 100;    % Save an image every N steps
num_spots = 10;         % Number of initial large spots
use_jet = true;         % true = jet (color), false = gray (grayscale)

% --- Uncertainty model for error propagation ---
sigma_F = 5e-4;         % Standard uncertainty of F
sigma_k = 5e-4;         % Standard uncertainty of k
sigma_num = 1e-6;       % Baseline numerical error floor
max_dt = 0.25;          % Internal substep upper bound for explicit Euler
clamp_max = 1.20;       % Physical concentration clamp to improve robustness

% --- List of parameters (F, k) ---
% Uncomment the lines you want to run and comment out the ones you don’t.
% The program will process all active lines and save results in separate folders
param_list = [
%%|   F     |   k    |   Pattern type
%%    0.010  0.050;  % isolated spots
%%    0.014  0.047;  % waves
%%    0.025  0.052;  % pulsating lotus
%%    0.022  0.059;  % cell war
      0.030  0.057;  % maze, Blaupunkt speaker pattern
%%    0.025  0.054;  % fingerprint 1
%%    0.038  0.060;  % fingerprint 2
%%    0.022  0.051;  % stripes and holes
%%    0.046  0.063;  % worms
%%    0.030  0.060;  % soliton
%%    0.035  0.065;  % spots
%%    0.037  0.057;  % firestorm
%%    0.095  0.056;  % soap bubbles
%%    0.082  0.059;  % puffball
%%    0.082  0.060;  % amoeboid
%%    0.062  0.061;  % frozen synergetics
%%    0.058  0.065;  % rods and loops
%%    0.042  0.059;  % concentric waves
];

try
    for idx = 1:rows(param_list)
        F = param_list(idx, 1);
        k = param_list(idx, 2);

        % --- Create folder for current parameter set ---
        foldername = sprintf("GrayScott_F%.3f_k%.3f", F, k);
        if ~exist(foldername, "dir")
            mkdir(foldername);
        end

        % --- Initial conditions ---
        U = ones(n, n);
        V = zeros(n, n);

        % Sensitivity states for linearized error propagation
        dU_dF = zeros(n, n);
        dV_dF = zeros(n, n);
        dU_dk = zeros(n, n);
        dV_dk = zeros(n, n);

        % History vectors for graphical post-analysis
        hist_step = [];
        hist_rms_sigmaV = [];
        hist_meanU = [];
        hist_meanV = [];

        % Old center disturbance initialization (optional):
##        r = 20; cx = n/2; cy = n/2;
##        U(cx-r:cx+r, cy-r:cy+r) = 0.50;
##        V(cx-r:cx+r, cy-r:cy+r) = 0.25 + 0.05*rand(r*2+1, r*2+1);

        % Generate random large spots
        for i = 1:num_spots
          radius = randi([5, 15]);               % Random radius
          cx = randi([radius+1, n-radius-1]);    % Center X
          cy = randi([radius+1, n-radius-1]);    % Center Y

          % Create circular mask
          [X, Y] = meshgrid(1:n, 1:n);
          mask = (X - cx).^2 + (Y - cy).^2 <= radius^2;

          % Reduce U and increase V inside the circle
          U(mask) = 0.50;
          V(mask) = 0.25 + 0.05 * rand();        % Small random variation
        end

        fprintf("\n▶ Running: F = %.3f, k = %.3f (%d/%d)\n", ...
             F, k, idx, rows(param_list));

        % --- Simulation loop ---
        for t = 0:steps
            [U, V, dU_dF, dV_dF, dU_dk, dV_dk, sigmaU_map, sigmaV_map] = ...
                grayscott_step(U, V, dU_dF, dV_dF, dU_dk, dV_dk, Du, Dv, ...
                F, k, dt, sigma_F, sigma_k, sigma_num, max_dt, clamp_max);

            hist_step(end+1) = t; %#ok<AGROW>
            hist_rms_sigmaV(end+1) = sqrt(mean(sigmaV_map(:).^2)); %#ok<AGROW>
            hist_meanU(end+1) = mean(U(:)); %#ok<AGROW>
            hist_meanV(end+1) = mean(V(:)); %#ok<AGROW>

            % Check for instability
            if any(isnan(U(:))) || any(isnan(V(:))) || ...
               any(isinf(U(:))) || any(isinf(V(:)))
                fprintf(["⚠ Instability detected at step %d! " ...
                    "Skipping to next parameter set.\n"], t);
                break; % Exit current simulation
            end

            % Visualization and saving
            if mod(t, save_interval) == 0
                if range(V(:)) < 1e-300  % Nearly homogeneous
                    fprintf(["⚠ V matrix is almost homogeneous, " ...
                        "skipping visualization.\n"]);
                else
                    if ~ishandle(1)
                        figure(1);
                    else
                        clf;
                    end

                    subplot(2,2,1);
                    imagesc(V);
                    if use_jet            % Select color palette
                        colormap(gca, jet);
                    else
                        colormap(gca, gray);
                    end
                    axis image; colorbar;
                    title(sprintf("V concentration (step: %d)", t), "fontsize", 12);
                    xlabel(sprintf("F = %.3f, k = %.3f", F, k), "fontsize", 10);

                    subplot(2,2,2);
                    imagesc(sigmaV_map);
                    colormap(gca, hot);
                    axis image; colorbar;
                    title("Propagated uncertainty of V", "fontsize", 12);
                    xlabel(sprintf("RMS(\sigma_V)=%.3e", hist_rms_sigmaV(end)), ...
                        "fontsize", 10);

                    subplot(2,2,3);
                    semilogy(hist_step, hist_rms_sigmaV, "LineWidth", 1.3);
                    grid on;
                    xlim([0, steps]);
                    xlabel("Step"); ylabel("RMS propagated error in V");
                    title("Error propagation over time", "fontsize", 12);

                    subplot(2,2,4);
                    plot(hist_step, hist_meanU, "b", "LineWidth", 1.2); hold on;
                    plot(hist_step, hist_meanV, "r", "LineWidth", 1.2); hold off;
                    grid on;
                    xlim([0, steps]);
                    xlabel("Step"); ylabel("Spatial mean concentration");
                    legend("mean(U)", "mean(V)", "Location", "best");
                    title("Global concentration trends", "fontsize", 12);

                    drawnow;

                    % Save a full analysis dashboard snapshot
                    filename = sprintf("%s/step_%05d.png", foldername, t);
                    print(filename, "-dpng", "-r120");
                    fprintf("  Saved: %s\n", filename);
                end
            end
        end

        % Save a compact result summary plot at the end for this parameter set
        summary_name = sprintf("%s/summary.png", foldername);
        figure(2); clf;
        subplot(1,2,1);
        semilogy(hist_step, hist_rms_sigmaV, "k", "LineWidth", 1.4); grid on;
        xlabel("Step"); ylabel("RMS propagated error in V");
        title("Error propagation summary");
        subplot(1,2,2);
        plot(hist_step, hist_meanU, "b", hist_step, hist_meanV, "r", "LineWidth", 1.2);
        grid on;
        xlabel("Step"); ylabel("Spatial mean concentration");
        legend("mean(U)", "mean(V)", "Location", "best");
        title("Mean concentration evolution");
        print(summary_name, "-dpng", "-r150");
        fprintf("  Saved: %s\n", summary_name);

        % --- Create GIF from PNG files ---
        try
            gif_name = sprintf("%s/animation.gif", foldername);
            img_files = dir(fullfile(foldername, "step_*.png"));
            if ~isempty(img_files)
                [~, order] = sort({img_files.name});
                img_files = img_files(order);
                fprintf("  ▶ Creating GIF: %s\n", gif_name);
                for i = 1:length(img_files)
                    img = imread(fullfile(foldername, img_files(i).name));
                    [A, map] = rgb2ind(img);
                    if i == 1
                        imwrite(A, map, gif_name, "gif", "LoopCount", Inf, ...
                            "DelayTime", 0.5);
                    else
                        imwrite(A, map, gif_name, "gif", "WriteMode", ...
                            "append", "DelayTime", 0.5);
                    end
                end
                fprintf("  GIF created: %s\n", gif_name);
            else
                fprintf("  No PNG files found, skipping GIF creation.\n");
            end
        catch gif_error
            fprintf("⚠ GIF creation error: %s\n", gif_error.message);
        end
    end
catch % Simulation can be interrupted anytime with Ctrl+C, then this block runs
    if ~ishandle(1)
        figure(1);
    else
        clf;
    end
    imagesc(V);
%%   colormap(jet);
    colormap(gray);
    axis image;
    colorbar;
    title("Simulation interrupted");
    xlabel(sprintf("F = %.3f, k = %.3f", F, k));
    drawnow;
    fprintf(["\n⚠ Simulation interrupted: parameter set %d (F = %.3f, " ...
      "k = %.3f), step: %d.\n"], idx, F, k, t);
    interrupted_file = sprintf("%s/interrupted_%05d.png", foldername, t);
    print(interrupted_file, "-dpng");
    fprintf("⚠ Current state saved: %s\n", interrupted_file);
end
