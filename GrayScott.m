%************************************************************************
% Gray–Scott model 2D simulation in GNU Octave
% Created by HA2ERZ & ChatGPT in 2025
% Released under CC0 1.0 Universal (Public Domain Dedication)
% You may copy, modify, distribute and perform the work,
% even for commercial purposes, without asking permission.
%************************************************************************

clear; clc;
pkg load image;

% --- Basic settings ---
n = 300;                % Grid size (tested up to 1200)
Du = 0.16;              % Diffusion coefficient for U
Dv = 0.08;              % Diffusion coefficient for V
dt = 1.0;               % Time step
steps = 5000;           % Number of iterations
save_interval = 100;    % Save an image every N steps
num_spots = 10;         % Number of initial large spots
use_jet = true;         % true = jet (color), false = gray (grayscale)

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

% --- Laplacian operator function ---
function L = laplacian(M)
    L = -4*M + circshift(M,[1,0]) + circshift(M,[-1,0]) ...
             + circshift(M,[0,1]) + circshift(M,[0,-1]);
end

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
            U_new = U + (Du*laplacian(U) - U.*V.^2 + F*(1 - U)) * dt;
            V_new = V + (Dv*laplacian(V) + U.*V.^2 - (F + k)*V) * dt;
            U = U_new; V = V_new;

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
                    imagesc(V);
                    if use_jet            % Select color palette
                      colormap(jet);
                    else
                      colormap(gray);
                    end
                    axis image;
                    colorbar;
                    title(sprintf("Gray-Scott model (step: %d)", t), ...
                       "fontsize", 14);
                    xlabel(sprintf("F = %.3f, k = %.3f", F, k), "fontsize", 12);
                    drawnow;
                    % Save PNG directly from matrix (1 cell = 1 pixel)
                    filename = sprintf("%s/step_%05d.png", foldername, t);
                    % All PNG files are saved in RGB format;
                    % only the applied colormap changes
                    if use_jet
                       cmap = jet(256);
                     else
                       cmap = gray(256);
                    end
                    img_rgb = ind2rgb(gray2ind(mat2gray(V), 256), cmap);
                    imwrite(img_rgb, filename);
                    % To save a screenshot, comment out the previous
                    % lines and uncomment the line below:
                    % print(filename, "-dpng");
                    fprintf("  Saved: %s\n", filename);
                end
            end
        end

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
