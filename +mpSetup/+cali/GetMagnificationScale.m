function [ transformations ] = GetMagnificationScale(cam1, cam2, cam3, cam4, inFocus1, inFocus2)
    SimilarityScore = [];
    
    for chIdx = 1:length(inFocus1)
        Plane9infocus = [];
        Plane1infocus = [];

        idxPlane1 = [inFocus1.globalch]==chIdx;
        idxPlane9 = [inFocus2.globalch]==chIdx;

        focus1     = inFocus1(idxPlane1).frame;
        focus9     = inFocus2(idxPlane9).frame;
    
        camIdx1 = inFocus1(idxPlane1).cam;
        camIdx9 = inFocus2(idxPlane9).cam;

        if camIdx1 == 1
            Plane1infocus = double(cam1(:,:,inFocus1(idxPlane1).ch,focus1));
        else
            Plane1infocus = double(cam2(:,:,inFocus1(idxPlane1).ch,focus1));
        end

        if camIdx9 == 3
            Plane9infocus = double(cam3(:,:,inFocus2(idxPlane9).ch,focus9));
        else
            Plane9infocus = double(cam4(:,:,inFocus2(idxPlane9).ch,focus9));
        end
    
        se = strel('disk', 10);
        bgPlane9 = imopen(Plane9infocus, se);
        bgPlane1 = imopen(Plane1infocus, se);
        Plane9infocus = Plane9infocus - bgPlane9;
        Plane1infocus = Plane1infocus - bgPlane1;

        figure()
        subplot(1,2,1)
        imshowpair(Plane1infocus,Plane9infocus)
        title("before correction");
        hold on

        config = "multimodal";
        transf = "similarity"; %similarity

        %% new approach
        scaleLimits = [0, 1.5];
        rotationLimits = [-1.2, 1.2]; % Rotation in degrees
        translationLimits = [-140, 140; 140, 140]; % Translation limits [X; Y]
        [optimizer, metric] = imregconfig("multimodal");
        optimizer.MaximumIterations  = 1000;
        tform_initial = imregtform(Plane9infocus, Plane1infocus, "similarity", optimizer, metric);
        T = tform_initial.T;
        scale_init = sqrt(T(1,1)^2 + T(2,1)^2); % Scale
        rotation_init = atan2d(T(2,1), T(1,1)); % Rotation (degrees)
        translation_init = T(3,1:2); % Translation (x, y)
        [best_tform, best_similarity] = customOptimizer(Plane9infocus, Plane1infocus, scale_init, rotation_init, translation_init, scaleLimits, rotationLimits, translationLimits, optimizer);
        transformations{chIdx, 1} = best_tform;
        transformations{chIdx, 2} = best_similarity;
        movingRegistered = imwarp(Plane9infocus, best_tform, "OutputView", imref2d(size(Plane1infocus)));

        % [optimizer,metric] = imregconfig(config);
        % optimizer.MaximumIterations = 10000;
        % tform = imregtform(Plane9infocus,Plane1infocus,transf, optimizer, metric);
        % transformations{chIdx, 1} = tform;
        % movingRegistered = imwarp(Plane9infocus,tform,"OutputView",imref2d(size(Plane1infocus)));
        % SimilarityScore(chIdx, 1) = multissim(movingRegistered,Plane1infocus);
        % transformations{chIdx, 2} = SimilarityScore(chIdx, 1);
        % transformations{chIdx, 3} = config;
        % transformations{chIdx, 4} = transf;

        subplot(1,2,2)
        imshowpair(Plane1infocus,movingRegistered);
        title("after correction");
        sgtitle(append("Plane ", num2str(chIdx), " x Plane ", num2str(chIdx+8)));
    end
end

% Custom Optimizer Function
function [best_tform, best_similarity] = customOptimizer(Plane9infocus, Plane1infocus, scale_init, rotation_init, translation_init, scaleLimits, rotationLimits, translationLimits, optimizer)
    best_similarity = -Inf; % Start with a very low similarity
    best_tform = []; % Will store the best transformation

    % Initial guess for the optimization parameters
    x0 = [scale_init, rotation_init, translation_init];

    % Get the optimization bounds
    lb = [scaleLimits(1), rotationLimits(1), translationLimits(1,1), translationLimits(2,1)];
    ub = [scaleLimits(2), rotationLimits(2), translationLimits(1,2), translationLimits(2,2)];

    % Define optimization options
    options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', optimizer.MaximumIterations);

    % Run the optimization loop for the defined maximum iterations
    for iteration = 1:optimizer.MaximumIterations
        % Run fmincon to optimize the transformation parameters
        x_opt = fmincon(@(x) costFunction(x, Plane9infocus, Plane1infocus), x0, [], [], [], [], lb, ub, [], options);

        % Extract optimized parameters
        scale = x_opt(1);
        theta = deg2rad(x_opt(2)); % Convert degrees to radians
        tx = x_opt(3);
        ty = x_opt(4);

        % Create transformation matrix
        T_restricted = [scale * cos(theta), -scale * sin(theta), 0;
                        scale * sin(theta), scale * cos(theta), 0;
                        tx, ty, 1];

        tform = affine2d(T_restricted);

        % Apply the transformation and compute the similarity score
        movingRegistered = imwarp(Plane9infocus, tform, "OutputView", imref2d(size(Plane1infocus)));
        similarity = multissim(movingRegistered, Plane1infocus); % Replace with your metric

        % If this transformation gives a better similarity, save it
        if similarity > best_similarity
            best_similarity = similarity;
            best_tform = tform;
        end

        % Update initial guess (for next iteration)
        x0 = x_opt; % Use the current solution as the starting point for next iteration
    end
end

% Cost Function to Minimize (based on similarity metric)
function error = costFunction(x, moving, fixed)
    scale = x(1);
    theta = deg2rad(x(2)); % Convert degrees to radians
    tx = x(3);
    ty = x(4);

    % Construct the transformation matrix
    T = [scale * cos(theta), -scale * sin(theta), 0;
         scale * sin(theta), scale * cos(theta), 0;
         tx, ty, 1];

    tform = affine2d(T);

    % Apply the transformation to the moving image
    movingTransformed = imwarp(moving, tform, "OutputView", imref2d(size(fixed)));

    % Compute the error (similarity score, lower error = better)
    error = immse(double(movingTransformed), double(fixed));
end