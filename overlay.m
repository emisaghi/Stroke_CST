% Overlay for CST in Stroke Project - Version June 9, 2016

clear all
close all
clc
dbstop if error

% Obtain file names for the maps, the tract files and the masks
FA_maps = dir('FA_maps');
FA_maps_name = {FA_maps(3:end).name};
AD_maps = dir('AD_maps');
AD_maps_name = {AD_maps(3:end).name};
ADC_maps = dir('ADC_maps');
ADC_maps_name = {ADC_maps(3:end).name};
RD_maps = dir('RD_maps');
RD_maps_name = {RD_maps(3:end).name};
Edema_masks = dir('Edema_masks');
Edema_masks_name = {Edema_masks(3:end).name};
% Hematoma_masks = dir('Hematoma_masks');
% Hematoma_masks_name = {Hematoma_masks(3:end).name};
Ipsi_tracts = dir('Ipsi_tracts');
Ipsi_tracts_name = {Ipsi_tracts(3:end).name};

% Make sure the files are stored in the same order and name in all
% subdirectories
FA_maps_only_name = regexprep(FA_maps_name,'.nii','');
AD_maps_only_name = regexprep(AD_maps_name,'.nii','');
ADC_maps_only_name = regexprep(ADC_maps_name,'.nii','');
RD_maps_only_name = regexprep(RD_maps_name,'.nii','');
Edema_masks_only_name = regexprep(Edema_masks_name, '.nii', '');
Ipsi_tracts_only_name = regexprep(Ipsi_tracts_name, '.mat', '');
% Hematoma_masks_only_name = regexprep(Hematoma_masks_name, '.nii', '');
if ~(isequal(FA_maps_only_name, Edema_masks_only_name))
   error('Please check FA maps and Edema masks file names!')
elseif ~(isequal(FA_maps_only_name, Ipsi_tracts_only_name))
    error('Please check FA maps and Ipsilateral tracts file names!')
elseif ~(isequal(FA_maps_only_name, AD_maps_only_name))
    error('Please check FA and AD maps file names!')
elseif ~(isequal(FA_maps_only_name, ADC_maps_only_name))
    error('Please check FA and ADC maps file names!')
elseif ~(isequal(FA_maps_only_name, RD_maps_only_name))
    error('Please check FA and RD maps file names!')
% elseif ~(isequal(FA_maps_only_name, Hematoma_masks_only_name))
%      error('Please check FA maps and Hematoma masks file names!')
else
    % Initialize the structure
    data(size(FA_maps_name, 2)).File_Name = '';
    data(size(FA_maps_name, 2)).Number_of_voxels = '';
    data(size(FA_maps_name, 2)).weighted_mean_FA = '';
    data(size(FA_maps_name, 2)).weighted_mean_AD = '';
    data(size(FA_maps_name, 2)).weighted_mean_ADC = '';
    data(size(FA_maps_name, 2)).weighted_mean_RD = '';
%     data(size(FA_maps_name, 2)).Hematoma_FA = '';
    data(size(FA_maps_name, 2)).Edema_FA = '';
    data(size(FA_maps_name, 2)).Edema_AD = '';
    data(size(FA_maps_name, 2)).Edema_ADC = '';
    data(size(FA_maps_name, 2)).Edema_RD = '';
    data(size(FA_maps_name, 2)).Min_Tract_Dist_From_Edema_Edges = '';
    data(size(FA_maps_name, 2)).count_1s = '';
    data(size(FA_maps_name, 2)).count_2s = '';
    data(size(FA_maps_name, 2)).count_3s = '';
    data(size(FA_maps_name, 2)).count_over_3s = '';
    data(size(FA_maps_name, 2)).FA_1s = '';
    data(size(FA_maps_name, 2)).FA_2s = '';
    data(size(FA_maps_name, 2)).FA_3s = '';
    data(size(FA_maps_name, 2)).FA_over_3s = '';

    for i = 1:size(FA_maps_name, 2)
        FA = E_DTI_read_nifti_file(strcat('FA_maps/', char(FA_maps_name(i)))); % Read FA map
        AD = E_DTI_read_nifti_file(strcat('AD_maps/', char(AD_maps_name(i)))); % Read AD map
        ADC = E_DTI_read_nifti_file(strcat('ADC_maps/', char(ADC_maps_name(i)))); % Read ADC map
        RD = E_DTI_read_nifti_file(strcat('RD_maps/', char(RD_maps_name(i)))); % Read RD map
        load(strcat('Ipsi_tracts/', char(Ipsi_tracts_name(i)))); % Read Tract file
        ROIMask = E_DTI_read_nifti_file(strcat('Edema_masks/', char(Edema_masks_name(i)))); % Read Edema mask
        ROIMask(ROIMask > 0) = 1;
        Tract_ROI_Overlay = ROIMask .* double(TractMask); % Calculate the number of streamlines passing through each of the voxels in the ROI
        FA_ROI = FA .* Tract_ROI_Overlay;
        mean_FA = sum(FA_ROI(:))/sum(Tract_ROI_Overlay(:)); % Calculate the weighted average of the FA values
        AD_ROI = AD .* Tract_ROI_Overlay;
        mean_AD = sum(AD_ROI(:))/sum(Tract_ROI_Overlay(:)); % Calculate the weighted average of the AD values
        ADC_ROI = ADC .* Tract_ROI_Overlay;
        mean_ADC = sum(ADC_ROI(:))/sum(Tract_ROI_Overlay(:)); % Calculate the weighted average of the ADC values
        RD_ROI = RD .* Tract_ROI_Overlay;
        mean_RD = sum(RD_ROI(:))/sum(Tract_ROI_Overlay(:)); % Calculate the weighted average of the RD values
        Tract_mask = double(TractMask > 0); % Build the binary tract mask
        Edema_FA = ROIMask .* Tract_mask .* FA; % Find the FA values of the intersection for Edema, Tract and the brain
        Edema_AD = ROIMask .* Tract_mask .* AD; % Find the AD values of the intersection for Edema, Tract and the brain
        Edema_ADC = ROIMask .* Tract_mask .* ADC; % Find the ADC values of the intersection for Edema, Tract and the brain
        Edema_RD = ROIMask .* Tract_mask .* RD; % Find the RD values of the intersection for Edema, Tract and the brain
        % Obtain nonzero values
        [~, ~, Edema_FA_nzv] = find(Edema_FA);
        [~, ~, Edema_AD_nzv] = find(Edema_AD);
        [~, ~, Edema_ADC_nzv] = find(Edema_ADC);
        [~, ~, Edema_RD_nzv] = find(Edema_RD);
        [d1, d2, d3]= ind2sub(size(ROIMask), find(ROIMask)); % Find the indices of the voxels that are inside the Edema mask
        d = [d1 d2 d3]; % Build a matrix of those indices
        [q1, q2, q3]= ind2sub(size(ROIMask .* Tract_mask), find(ROIMask .* Tract_mask)); % Find the indices of the intersection voxels
        q = [q1 q2 q3]; % Build a matrix of those indices
        o = zeros(size(q1, 1), 6);
        for j = 1:size(q1) % Find the direct distance from the edges of the mask
            r1 = d(d2 == q2(j) & d3 == q3(j), 1);
            r2 = d(d1 == q1(j) & d3 == q3(j), 2);
            r3 = d(d1 == q1(j) & d2 == q2(j), 3);
            o(j, :) = [1+min(q1(j)-min(r1), max(r1)-q1(j)), 1+min(q2(j)-min(r2), max(r2)-q2(j)), 1+min(q3(j)-min(r3), max(r3)-q3(j)), q1(j), q2(j), q3(j)];
        end
        res = [min(o(:, 1:3), [], 2), o(:, 4:6)];
        %indices
        a1 = res(res(:, 1) == 1, 2:4);
        a2 = res(res(:, 1) == 2, 2:4);
        a3 = res(res(:, 1) == 3, 2:4);
        a4 = res(res(:, 1) > 3, 2:4);
        FA1 = zeros(size(a1, 1), 1);
        FA2 = zeros(size(a2, 1), 1);
        FA3 = zeros(size(a3, 1), 1);
        FA4 = zeros(size(a4, 1), 1);
        for k = 1:size(a1, 1)
            FA1(k) = FA(a1(k, 1), a1(k, 2), a1(k, 3));
        end
        for k = 1:size(a2, 1)
            FA2(k) = FA(a2(k, 1), a2(k, 2), a2(k, 3));
        end
        for k = 1:size(a3, 1)
            FA3(k) = FA(a3(k, 1), a3(k, 2), a3(k, 3));
        end
        for k = 1:size(a4, 1)
            FA4(k) = FA(a4(k, 1), a4(k, 2), a4(k, 3));
        end
        % Count the number of voxels that have a distance of 1, 2, 3 or
        % greater than 3 number of voxels from the edge of the mask
        count1 = nnz(res(:, 1) == 1);
        count2 = nnz(res(:, 1) == 2);
        count3 = nnz(res(:, 1) == 3);
        count4 = nnz(res(:, 1) > 3);
        
%         ROIMask = E_DTI_read_nifti_file(strcat('Hematoma_masks/', char(Hematoma_masks_name(i)))); % Read Hematoma mask
%         ROIMask(ROIMask > 0) = 1; % The hematoma mask is composed of 0s and 2s. Change 2s to 1s.
%         Hematoma_FA = ROIMask .* Tract_mask .* FA; % Find the intersection for Hematoma, Tract and the brain
%         [~, ~, Hematoma_FA_nzv] = find(Hematoma_FA); % Obtain nonzero values
%         data(i) = struct('File_Name', char(FA_maps_name(i)), 'Hematoma_FA', Hematoma_FA_nzv, 'Edema_FA', Edema_FA_nzv); % Fill the structure with the obtained values
        data(i) = struct('File_Name', char(FA_maps_name(i)), 'Number_of_voxels', size(Edema_FA_nzv, 1), 'weighted_mean_FA', mean_FA, 'weighted_mean_AD', mean_AD, 'weighted_mean_ADC', mean_ADC, 'weighted_mean_RD', mean_RD, 'Edema_FA', Edema_FA_nzv, 'Edema_AD', Edema_AD_nzv, 'Edema_ADC', Edema_ADC_nzv, 'Edema_RD', Edema_RD_nzv, 'Min_Tract_Dist_From_Edema_Edges', res, 'count_1s', count1, 'count_2s', count2, 'count_3s', count3, 'count_over_3s', count4, 'FA_1s', FA1, 'FA_2s', FA2, 'FA_3s', FA3, 'FA_over_3s', FA4); % Fill the structure with the obtained values
    end

    % Save the results into a spreadsheet
    save('data.mat', 'data');
    t = struct2table(data);
    writetable(t, 'data.txt');
    writetable(t, 'data.xlsx');

end