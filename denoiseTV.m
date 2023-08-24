function u = denoiseTV(f, lambda, numIterations)

    % Convert input image to double

    f = double(f);

    

    % Normalize the image to the range [0, 1]

    f = (f - min(f(:))) / (max(f(:)) - min(f(:)));

    

    % Initialize the denoised image as the input image

    u = f;

    

    % Compute gradient operators

    dx = [-1 1];

    dy = [-1; 1];

    

    % Perform Total Variation denoising

    for iter = 1:numIterations

        % Compute gradients

        ux = imfilter(u, dx, 'replicate');

        uy = imfilter(u, dy, 'replicate');

        

        % Compute divergence

        div = divergence(ux, uy);

        

        % Update the denoised image

        u = u + lambda * div;

    end

    

    % Rescale the denoised image back to the original intensity range

    u = u * (max(f(:)) - min(f(:))) + min(f(:));

end

