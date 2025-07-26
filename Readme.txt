===========================================================
Gray–Scott Reaction-Diffusion Simulation (GNU Octave)
===========================================================
Created by HA2ERZ & ChatGPT in 2025
-----------------------------------------------------------
Overview:
-----------------------------------------------------------
This script simulates the Gray–Scott reaction-diffusion system
in 2D, producing complex self-organizing patterns such as
spots, stripes, labyrinths, worms, and more.

Each selected parameter set creates:
 - Multiple PNG snapshots (grid-based, 1 cell = 1 pixel)
 - An animated GIF showing pattern evolution

-----------------------------------------------------------
Features:
-----------------------------------------------------------
- Multiple (F, k) parameter sets can run automatically
- Outputs exact-scale PNG images (not distorted)
- RGB PNG saving with selectable colormap:
    - Jet (colorful)
    - Gray (grayscale)
- Large random spots initialization for faster pattern growth
- Instability detection: skips unstable runs
- Automatic GIF animation from generated PNG frames

-----------------------------------------------------------
Requirements:
-----------------------------------------------------------
- GNU Octave 7.x or later
- Image package:
    pkg load image

-----------------------------------------------------------
How to Use:
-----------------------------------------------------------
1. Open GrayScott.m in GNU Octave.

2. Set basic parameters at the beginning of the script:
    n = 500;            % Grid size
    steps = 15000;      % Total simulation steps
    save_interval = 1000; % Save every N steps
    num_spots = 50;     % Number of initial random spots

3. Choose colormap mode:
    use_jet = true;    % true = jet (color), false = gray (grayscale)

4. Select parameter sets by editing `param_list`:
    param_list = [
        0.035 0.065;    % spots
        0.030 0.057;    % labyrinth
    ];
   - Enable by uncommenting
   - Disable by commenting out

5. Run the script:
    octave GrayScott.m

-----------------------------------------------------------
Output:
-----------------------------------------------------------
Each parameter set will create a separate folder:
    GrayScott_F0.035_k0.065/
        step_00000.png
        step_01000.png
        ...
        animation.gif

PNG files:
    - Are always RGB (3-channel) images.
    - The applied colormap (jet or gray) changes only the colors.

GIF file:
    - Created from the PNG frames automatically.
    - Matches the selected colormap.

-----------------------------------------------------------
Visualization:
-----------------------------------------------------------
- Preview during simulation:
    imagesc(V);
    if use_jet
        colormap(jet);
    else
        colormap(gray);
    end
    axis image;
    colorbar;

- PNG saving in the script:
    if use_jet
        cmap = jet(256);
    else
        cmap = gray(256);
    end
    img_rgb = ind2rgb(gray2ind(mat2gray(V), 256), cmap);
    imwrite(img_rgb, filename);

-----------------------------------------------------------
Initialization:
-----------------------------------------------------------
Default: Large random spots of different sizes.
(Old central disturbance method is commented out in the code.)

To change:
 - Modify `num_spots` for density.
 - Change radius range in initialization for larger or smaller spots.

-----------------------------------------------------------
Example Parameter Sets:
-----------------------------------------------------------
| Pattern        | F       | k       |
|--------------- |---------|---------|
| Isolated spots | 0.010   | 0.050   |
| Waves          | 0.014   | 0.047   |
| Labyrinth      | 0.030   | 0.057   |
| Spots          | 0.035   | 0.065   |
| Worms & Loops  | 0.058   | 0.065   |
| Soap Bubbles   | 0.095   | 0.056   |

-----------------------------------------------------------
License:
-----------------------------------------------------------
CC0 1.0 Universal (Public Domain Dedication)
You are free to copy, modify, distribute, and use for any purpose,
even commercially, without asking permission.
