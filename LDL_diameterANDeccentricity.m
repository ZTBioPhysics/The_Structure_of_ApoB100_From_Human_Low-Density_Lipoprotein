% Load image stack
stack = ReadMRC('LDL_200Classes.mrc');
num_images = length(stack(1,1,:));
apix = 450/(length(stack(:,1,1)))*1.09;

% Initialize variables for storing diameter and eccentricity
diameter = zeros(1,num_images);
average_diameter = zeros(1,num_images);
eccentricity = zeros(1,num_images);
major_lengths = zeros(1,num_images);
minor_lengths = zeros(1,num_images);

% Loop through images
for i = 1:num_images
    
    % Load image and convert to binary mask
    image = stack(:,:,i);
    bw = imbinarize(image);
    
    % Find connected components in mask
    cc = bwconncomp(bw);
    
    % Find largest connected component
    num_pixels = cellfun(@numel, cc.PixelIdxList);
    [~, max_idx] = max(num_pixels);
    bw(:) = false;
    bw(cc.PixelIdxList{max_idx}) = true;
    
    % Get region properties of connected component
    props = regionprops(bw, 'MajorAxisLength','MinorAxisLength','Eccentricity','Centroid','Orientation');
    diameter(i) = props.MajorAxisLength*apix;
    eccentricity(i) = props.Eccentricity;
    major_lengths(i) = props.MajorAxisLength*apix;
    minor_lengths(i) = props.MinorAxisLength*apix;
    average_diameter(i) = (props.MajorAxisLength+props.MinorAxisLength)/2*apix;
    
    % Display image and add axis lines
    figure(1);
    subplot(1,2,1)
    imshow(image);
    subplot(1,2,2)
    imshow(bw)
    hold on;
    plot(props.Centroid(1), props.Centroid(2), 'k+', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Draw major axis line in red
    x1 = props.Centroid(1) + 0.5*props.MajorAxisLength*cosd(-props.Orientation);
    y1 = props.Centroid(2) + 0.5*props.MajorAxisLength*sind(-props.Orientation);
    x2 = props.Centroid(1) - 0.5*props.MajorAxisLength*cosd(-props.Orientation);
    y2 = props.Centroid(2) - 0.5*props.MajorAxisLength*sind(-props.Orientation);
    line([x1 x2], [y1 y2], 'Color', 'r', 'LineWidth', 2);
    
    % Draw minor axis line in blue
    x1 = props.Centroid(1) + 0.5*props.MinorAxisLength*sind(-props.Orientation);
    y1 = props.Centroid(2) - 0.5*props.MinorAxisLength*cosd(-props.Orientation);
    x2 = props.Centroid(1) - 0.5*props.MinorAxisLength*sind(-props.Orientation);
    y2 = props.Centroid(2) + 0.5*props.MinorAxisLength*cosd(-props.Orientation);
    line([x1 x2], [y1 y2], 'Color', 'b', 'LineWidth', 2);
    
    % Add title and labels
    title(sprintf('Image %d: Diameter = %.2f, Eccentricity = %.2f', i, diameter, props.Eccentricity));
    
    % add a pause for visualization
    pause(0.05)

end

% histograms
figure(2);
subplot(1,2,1);
histogram(average_diameter);
xlabel('Average Diameter');
ylabel('Count');
title('Histogram of Diameter');

subplot(1,2,2);
histogram(eccentricity);
xlabel('Eccentricity');
ylabel('Count');
title('Histogram of Eccentricity');

% scatter plots
figure(3);
subplot(1,2,1);
scatter(major_lengths, minor_lengths);
xlabel('Major Axis Length');
ylabel('Minor Axis Length');

subplot(1,2,2);
scatter(average_diameter, eccentricity);
xlabel('Average Diameter');
ylabel('Eccentricity');
hold off;

%% cluster images with k-means based on major and minor axis lengths

% Perform k-means clustering
num_clusters = 10;
data = [major_lengths', minor_lengths'];
[idx, centroids] = kmeans(data, num_clusters);

% Plot scatter plot of Major and Minor axis lengths with cluster colors
figure;
gscatter(data(:,1), data(:,2), idx);
hold on;
scatter(centroids(:,1), centroids(:,2),100, 'kx', 'LineWidth', 2);
xlabel('Major Axis Length');
ylabel('Minor Axis Length');
title('K-means Clustering of Major and Minor Axis Lengths');
hold off;

% Find the representative image index for each cluster
rep_image_idx = zeros(1, num_clusters);
for i = 1:num_clusters
    [~, rep_image_idx(i)] = min(vecnorm(data - centroids(i,:), 2, 2));
end

% Display one representative image from each cluster
figure;
for i = 1:num_clusters
    subplot(ceil(sqrt(num_clusters)), ceil(sqrt(num_clusters)), i);
    rep_image = stack(:,:,rep_image_idx(i));
    imshow(rep_image);
    title(sprintf('Cluster %d', i));
end

% Print a list of all the images in each group
for i = 1:num_clusters
    group_images = find(idx == i);
    fprintf('Group %d: %s\n', i, mat2str(group_images));
end


%% group images baesd on their major axis lengths

% Divide the images into 10 groups based on their major axis lengths
num_groups = 10;
[~, ~, group_idx] = histcounts(major_lengths, num_groups);

% Create a custom colormap with 10 colors
colormap = lines(num_groups);

% Plot the histogram of major axis lengths with different color bins for the groups
figure;
hold on;
for i = 1:num_groups
    group_major_lengths = major_lengths(group_idx == i);
    histogram(group_major_lengths, 'FaceColor', colormap(i,:), 'BinMethod', 'integers', 'DisplayStyle', 'bar', 'EdgeColor', 'none');
end
xlabel('Major Axis Length');
ylabel('Count');
title('Histogram of Major Axis Lengths with Group Colors');
hold off;

% Display one representative image from each group
figure;
for i = 1:num_groups
    subplot(ceil(sqrt(num_groups)), ceil(sqrt(num_groups)), i);
    group_images = find(group_idx == i);
    rep_image = stack(:,:,group_images(1));
    imshow(rep_image);
    title(sprintf('Group %d', i));
end

% Print a list of all the images in each group
for i = 1:num_groups
    group_images = find(group_idx == i);
    fprintf('Group %d: %s\n', i, mat2str(group_images));
end


