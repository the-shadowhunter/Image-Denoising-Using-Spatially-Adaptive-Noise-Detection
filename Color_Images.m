% Main Script for Colored Images
I = imread("F:\Image_Denoising_Dataset\RGB_Images\Noisy_Images\20.png"); %Path_to_Noisy_Image
% I = imresize(I,[512,512]);
I = double(I);  % Convert to double for calculations
[rows, cols] = size(I);
% Check if the image is RGB
if size(I, 3) == 3
    % Split the RGB image into channels
    red_channel = I(:, :, 1);
    green_channel = I(:, :, 2);
    blue_channel = I(:, :, 3);

    % Process each channel
    noise_map_red = process_image(red_channel);
    noise_map_green = process_image(green_channel);
    noise_map_blue = process_image(blue_channel);

    denoised_red = restore_image(red_channel, noise_map_red, size(red_channel, 1), size(red_channel, 2));
    denoised_green = restore_image(green_channel, noise_map_green, size(green_channel, 1), size(green_channel, 2));
    denoised_blue = restore_image(blue_channel, noise_map_blue, size(blue_channel, 1), size(blue_channel, 2));

    denoised_red = postprocess_image(red_channel, denoised_red, noise_map_red);
    denoised_green = postprocess_image(green_channel, denoised_green, noise_map_green);
    denoised_blue = postprocess_image(blue_channel, denoised_blue, noise_map_blue);

    % Combine the processed channels into a single RGB image
    denoised_image = cat(3, denoised_red, denoised_green, denoised_blue);

    % Display the results
    figure;
    subplot(1, 2, 1);
    imshow(uint8(I));  % Original Image
    title('Original Image');
    subplot(1, 2, 2);
    imshow(uint8(denoised_image));  % Denoised Image
    title('Denoised Image');
else
    error('The input image is not an RGB image.');
end


%%
function [noise_map]= process_image(I)
    % I = imread(image_path);
    % if size(I, 3) == 3
    %     I = rgb2gray(I);  % Convert to grayscale if image is RGB
    % end
    I = double(I);  % Convert to double for calculation
    [rows, cols] = size(I);

    % Parameters
    noise_threshold = 30;
    intensity_threshold = 20;
    S_ref = 20;
    extreme_values = [0,255];

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
            max_window_size = min(min(min(i-1, rows-i), min(j-1, cols-j)),floor(sqrt(min(rows,cols))));
            % max_window_size = floor(sqrt(max(rows,cols)));

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
                    % noisy_pixel_found = true;
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
            max_window_size = max(min(min(i-1, rows-i), min(j-1, cols-j)),floor(sqrt(min(rows,cols))));

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
    

    % Post-processing step: Handle black pixels and refine the image
    for i = 2:rows_padded-1
        for j = 2:cols_padded-1
            if padded_image(i, j) == 0  % Check for black (noisy) pixels
                window_size = 1;
                valid_neighbors_found = false;
                max_postproc_window = max(min(min(i-1, rows_padded-i), min(j-1, cols_padded-j)),floor(sqrt(min(rows_padded,cols_padded))));% Maximum window size for post-processing

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
    denoised_image = medfilt2(denoised_image, [3 3]);

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

