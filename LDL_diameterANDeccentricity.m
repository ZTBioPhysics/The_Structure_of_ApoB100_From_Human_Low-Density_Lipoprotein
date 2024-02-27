% Load image stack
stack = ReadMRC('cryosparc_P5_J294_templates_selected.mrc');
num_images = length(stack(1,1,:));
apix = 450/(length(stack(:,1,1)))*1.09;

% Initialize variables for storing diameter and circularity
diameters = zeros(num_images, 1);
eccentricities = zeros(num_images, 1);

% Loop through images and compute diameter and circularity
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
    props = regionprops(bw, 'MajorAxisLength', 'MinorAxisLength', 'Area','Eccentricity');
    diameter = props.MajorAxisLength*apix;
    diameters(i) = diameter;
%     circularities(i) = 4*pi*props.Area/(props.MajorAxisLength^2);
    eccentricities(i) = props.Eccentricity;
end

% Compute pairwise distances between the images based on diameter and circularity
features = [diameters, eccentricities];
distances = pdist(features);

% Plot the pairwise distance matrix
figure;
imagesc(squareform(distances));
colorbar;
title('Pairwise Distance Matrix');

% Perform hierarchical clustering and plot the dendrogram
tree = linkage(features, 'average', 'euclidean');
figure;
dendrogram(tree, num_images);
title('Dendrogram');

% Cluster images using hierarchical clustering
num_clusters = 10;
cluster_indices = cluster(tree, 'maxclust', num_clusters);

% Compute the mean diameter and circularity for each cluster
D = zeros(num_clusters, 1);
E = zeros(num_clusters, 1);
for i = 1:num_clusters
    index = find(cluster_indices == i);
    if ~isempty(index)
        D(i) = mean(diameters(index));
        E(i) = mean(eccentricities(index));
    end
end

% Plot a representative image from each cluster
figure;
for i = 1:num_clusters
    % Find the index of the images in the current cluster
    index = find(cluster_indices == i);
    
    if isempty(index)
        % If the cluster is empty, skip to the next cluster
        continue
    end
    
    % Find the index of the image with the closest diameter to the mean diameter of the cluster
    [~, representative_index] = min(abs(D(i) - diameters(index)));
    representative_image = stack(:, :, index(representative_index));
    
    % Display the representative image in a subplot
    subplot(2, 5, i);
    imshow(representative_image);
    title(sprintf('Cluster %d', i));
    
    % Display the mean diameter and circularity of the cluster
    xlabel(sprintf('Diameter = %.2f\nEccentricity = %.2f', D(i), E(i)));
end



