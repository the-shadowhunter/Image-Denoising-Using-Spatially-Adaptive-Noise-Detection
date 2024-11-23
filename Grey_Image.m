path = "F:\Kaggle_2\Noisy_folder\noisy_45031049054_571119694c_c.jpg"; %Path_to_Noisy_Image
I = imread(path);
I = double(I);  % Convert to double for calculation
[rows, cols] = size(I);
map_noisy = process_image(path);
denoised_image = restore_image(I, map_noisy, rows, cols);
denoised_image = postprocess_image(I, denoised_image, map_noisy);
%%
% Specify the folder containing noisy images
sourceFolder = "F:\Image_Denoising_Dataset\Grey_Scale_Images\Noisy";
% Specify the folder to save the denoised images
destinationFolder = "F:\Image_Denoising_Dataset\Grey_Scale_Images\Denoised_Images";

% Create the destination folder if it does not exist
if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end

% Get a list of all image files in the source folder
imageFiles = dir(fullfile(sourceFolder, '*.*'));
validExtensions = {'.jpg', '.jpeg', '.png', '.bmp', '.tiff'}; % Add more formats if needed

% Loop through each file in the source folder
for i = 1:length(imageFiles)
    [~, name, ext] = fileparts(imageFiles(i).name); % Get the file name and extension
    if ismember(lower(ext), validExtensions) % Check if the file is a valid image
        % Read the noisy image
        noisyImagePath = fullfile(sourceFolder, imageFiles(i).name);
        noisyImage = imread(noisyImagePath);
        noisyImage = double(noisyImage); % Convert to double for processing
        
        % Get image dimensions
        [rows, cols, ~] = size(noisyImage);
        
        % Denoising process
        map_noisy = process_image(noisyImagePath); % Preprocess the image (custom function)
        denoisedImage = restore_image(noisyImage, map_noisy, rows, cols); % Restore the image (custom function)
        denoisedImage = postprocess_image(noisyImage, denoisedImage, map_noisy); % Postprocess (custom function)
        
        % Convert denoised image back to uint8 for saving
        denoisedImage = uint8(denoisedImage);
        
        % Save the denoised image in the destination folder with the same name
        denoisedImagePath = fullfile(destinationFolder, [name, ext]);
        imwrite(denoisedImage, denoisedImagePath);
        fprintf('Processed and saved: %s\n', denoisedImagePath);
    end
end

disp('All noisy images have been processed and saved as denoised images.');

%%
function [noise_map]= process_image(image_path)
    I = imread(image_path);
    % if size(I, 3) == 3
    %     I = rgb2gray(I);  % Convert to grayscale if image is RGB
    % end
    I = double(I);  % Convert to double for calculation
    [rows, cols] = size(I);

    % Parameters
    noise_threshold = 30;
    intensity_threshold = 20;
    S_ref = 20;
    extreme_values = [0, 255];

    % Initialize output image and noise map
    denoised_image = I;
    noise_map = zeros(rows, cols);  % 1 if noisy, 0 if non-noisy

    % Corrected function to check ambivalence
    is_ambivalent = @(pixel, neighbors_non_extreme, intensity_threshold) ...
        ismember(pixel, extreme_values) || all(abs(pixel - neighbors_non_extreme) > intensity_threshold);

    % Main loop for processing each pixel
    for i = 1:rows
        for j = 1:cols
            window_size = 1;
            noisy_pixel_found = false;
            max_window_size = min(min(min(i-1, rows-i), min(j-1, cols-j)), floor(sqrt(min(rows,cols))-3));

            while ~noisy_pixel_found && window_size <= max_window_size
                % Define window boundaries
                x_min = max(1, i - window_size);
                x_max = min(rows, i + window_size);
                y_min = max(1, j - window_size);
                y_max = min(cols, j + window_size);
                window = I(x_min:x_max, y_min:y_max);
                Ia = I(i, j);
                neighbors = window(:);  % Keep original neighbors (including extreme values)
                neighbors_non_extremes = neighbors(~ismember(neighbors, extreme_values));  % Exclude extremes
                if isempty(neighbors_non_extremes)
                    neighbors_non_extremes = neighbors;  % Fallback to original neighbors
                end
                current_pixel = Ia;

                % Step: Check for ambivalence using non-extreme neighbors
                ambivalent = is_ambivalent(Ia, neighbors_non_extremes, intensity_threshold);
                neighbors_1 = neighbors_non_extremes(neighbors_non_extremes ~= current_pixel);
                num_ambivalent = sum(arrayfun(@(p) is_ambivalent(p, neighbors_1, intensity_threshold), neighbors_non_extremes));
                num_pixels = numel(neighbors);
                app = (num_ambivalent / num_pixels) * 100;  % Ambivalent Pixel Percentage (APP)

                % Step: Noise detection based on APP
                if num_ambivalent == 1
                    noisy_pixel_found = true;
                    noise_map(i, j) = 1;
                    break;
                end

                if num_ambivalent > 1
                    % Calculate MSI if APP is insufficient for noise detection
                    wp_ref = round(0.3 * numel(neighbors));
                    SIR = abs(neighbors_non_extremes - Ia);
                    SIR_sorted = sort(SIR, 'ascend');
                    MSI = sum(SIR_sorted / S_ref) - wp_ref;

                    if app <= noise_threshold
                        if ismember(Ia, extreme_values)
                            noisy_pixel_found = true;
                            noise_map(i, j) = 1;
                        else
                            if MSI > 0
                                noisy_pixel_found = true;
                                noise_map(i, j) = 1;
                            else
                                if window_size >= max_window_size
                                    noisy_pixel_found = true;
                                    noise_map(i, j) = 0;
                                else
                                    window_size = window_size + 1;
                                end
                            end
                        end
                    else
                        if MSI > 0
                            noisy_pixel_found = true;
                            noise_map(i, j) = 1;
                        else
                            if window_size >= max_window_size
                                noisy_pixel_found = true;
                                noise_map(i, j) = 0;
                            else
                                window_size = window_size + 1;
                            end
                        end
                    end
                else
                    noise_map(i, j) = 0;
                    noisy_pixel_found = true;
                end
            end
        end
    end

    % Display original and denoised images
    figure;
    subplot(1, 2, 1);
    imshow(uint8(I));  % Original Image
    title('Original Image');
    subplot(1, 2, 2);
    imshow(noise_map, []);
    title('Noise Map');
end

%%
function [denoised_image] =  restore_image(I, noise_map, rows, cols)
    % Initialize variables
    map_noise = noise_map;
    denoised_image = I;
    is_noisy = @(pixel) ismember(pixel, 1);
    extreme_values = 1;

    % Loop through each pixel
    for i = 1:rows
        for j = 1:cols
            window_size = 1;
            noisy_pixel_restored = false;
            max_window_size = min(min(min(i-1, rows-i), min(j-1, cols-j)), floor(sqrt(min(rows,cols))-3));

            while ~noisy_pixel_restored && window_size <= max_window_size
                % Define window boundaries
                x_min = max(1, i - window_size);
                x_max = min(rows, i + window_size);
                y_min = max(1, j - window_size);
                y_max = min(cols, j + window_size);

                window_1 = I(x_min:x_max, y_min:y_max);
                window_2 = map_noise(x_min:x_max, y_min:y_max);
                Ia_1 = I(i, j);
                Ia_2 = map_noise(i, j);

                neighbors = window_1(:);
                neighbors_nature = window_2(:);
                neighbors_non_noisy = neighbors(~ismember(neighbors_nature, extreme_values));  % Exclude noisy pixels

                noisy = is_noisy(Ia_2);
                num_noisy = numel(neighbors) - numel(neighbors_non_noisy);
                num_pixels = numel(neighbors);
                np = (num_noisy / num_pixels) * 100;  % Noisy Pixel Percentage

                % Step: Noise detection and restoration
                if num_noisy == 1 || np < 30 || window_size >= max_window_size
                    denoised_image(i, j) = median(neighbors_non_noisy);
                    noisy_pixel_restored = true;
                else
                    window_size = window_size + 1;
                end
            end
        end
    end

    % Display original and denoised images
    figure;
    subplot(1, 2, 1);
    imshow(uint8(I));  % Original Image
    title('Original Image');
    subplot(1, 2, 2);
    imshow(uint8(denoised_image));  % Denoised Image
    title('Denoised Image');
end
%%
function denoised_image = postprocess_image(I, denoised_image, noise_map)
    % Mirror pad the original image to handle boundaries
    padded_image = padarray(denoised_image, [1, 1], 'symmetric'); % Padding with symmetric mirroring
    [rows_padded, cols_padded] = size(padded_image);
    [rows, cols] = size(I);

    % Post-processing step: Handle black pixels and refine the image
    max_postproc_window = floor(sqrt(min(rows,cols))-3);  % Maximum window size for post-processing
    for i = 2:rows_padded-1
        for j = 2:cols_padded-1
            if padded_image(i, j) == 0  % Check for black (noisy) pixels
                window_size = 1;
                valid_neighbors_found = false;

                while ~valid_neighbors_found && window_size <= max_postproc_window
                    % Define window boundaries in padded coordinates
                    x_min = max(1, i - window_size);
                    x_max = min(rows_padded, i + window_size);
                    y_min = max(1, j - window_size);
                    y_max = min(cols_padded, j + window_size);

                    % Extract neighbors
                    window = padded_image(x_min:x_max, y_min:y_max);
                    valid_neighbors = window(window > 0);  % Exclude black (noisy) pixels

                    if ~isempty(valid_neighbors)
                        % Replace black pixel with median of valid neighbors
                        padded_image(i, j) = median(valid_neighbors);
                        valid_neighbors_found = true;
                    else
                        % Increase window size if no valid neighbors are found
                        window_size = window_size + 1;
                    end
                end

                % Fallback: If no valid neighbors are found within max window size
                if ~valid_neighbors_found
                    padded_image(i, j) = I(i-1, j-1);  % Revert to original value (adjust for padding)
                end
            end
        end
    end

    % Remove padding to restore original image size
    denoised_image = padded_image(2:end-1, 2:end-1);

    % Optional: Apply a median filter for further smoothing
    % denoised_image = medfilt2(denoised_image, [3 3]);

    % Display results
    figure;
    subplot(1, 3, 1);
    imshow(uint8(I));  % Original Image
    title('Original Image');
    subplot(1, 3, 2);
    imshow(noise_map, []);
    title('Noise Map');
    subplot(1, 3, 3);
    imshow(uint8(denoised_image));  % Denoised Image
    title('Final Denoised Image');
end
