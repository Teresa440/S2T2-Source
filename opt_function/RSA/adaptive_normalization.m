function norm_arch_objs = adaptive_normalization(archive,nobjs)
%ADAPTIVE_NORMALIZATION: Normalizes points in the objective space w.r.t. the
%ideal point z_ideal
%   Normalization is performed using as extreme points to construct the
%   hyperplane the points that minimizes the Achievement Scalarizing 
%   Function (ASF), or (in case of duplicates or negative intercepts) the 
%   points corresponding to a max value in the objective space (in total 
%   they are nobjs).

points_id = archive.slot_status ~= 0;

z_ideal = min(archive.objs(points_id,:)); % ideal point
norm_archive_obj = archive.objs;
norm_archive_obj(points_id,:) = norm_archive_obj(points_id,:) - z_ideal; % offsetting the archive by ideal point

d = tril(ones(nobjs)*1e-6, -1) + triu(ones(nobjs)*1e-6, 1) + diag(ones(nobjs,1)); % identity matrix, with 1e-6 in place of the zeros
extreme_points = zeros(nobjs,nobjs);
f = norm_archive_obj(points_id,:);
for k = 1:nobjs
    dk = d(k,:);
    [~, id] = min(max(f./dk,[],2),[],1); % Achievement Scalarizing Function (ASF)
    extreme_points(k,:) = f(id,:); % points which minimizes ASF
    norm_archive_obj(points_id,k) = norm_archive_obj(points_id,k)/extreme_points(k,k);
end


[~,IA,~] = unique(extreme_points,'rows'); % check if there are duplicate extreme points 
if length(IA) < nobjs
    % if this is the case the procedure without ASF is followed
    disp("Invalid hyperplane using Achievement Scalarizing Function (duplicate extreme points)");
    norm_arch_objs = extreme_points_without_ASF(nobjs, norm_archive_obj, archive);
else
    plane_coeff = extreme_points\ones(nobjs,1); % determination of hyperplane coefficients
    intercepts = (1./plane_coeff)'; % calculating intercepts for every axis, transposed to make it a row vector    
    norm_arch_objs(points_id,:) = f./intercepts; % normalization of the archive
    % check if there are negative intercepts
    if any(intercepts < 0)
        % if this is the case the results are discarded and the procedure without ASF is followed
        disp("Invalid hyperplane using Achievement Scalarizing Function (negative intercepts)");
        norm_arch_objs = extreme_points_without_ASF(nobjs, norm_archive_obj, archive);
    end
end

end

%% Local Functions

function norm_arch_objs = extreme_points_without_ASF(nobjs, norm_archive_obj, archive)
%EXTREME_POINTS_WITHOUT_ASF: selects extreme points for hyperplane
%contruction
%   selects extreme points by chosing the ones which have max objective
%   value, without using ASF

extreme_points = zeros(nobjs,nobjs); % extreme point, used to construct the hyperplane
[~, id] = max(norm_archive_obj);
for i = 1:nobjs
    extreme_points(i,:) = archive.objs(id(i),:);
end

plane_coeff = extreme_points\ones(nobjs,1); % determination of hyperplane coefficients
intercepts = (1./plane_coeff)'; % calculating intercepts for every axis, transposed to make it a row vector
norm_arch_objs = archive.objs./intercepts; % normalization of the archive

end