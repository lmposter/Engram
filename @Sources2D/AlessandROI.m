function ind_del = AlessandROI(obj, vid_path, ind, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.
vid_path1 = fullfile(vid_path, 'msvideo.avi');
vid = VideoReader(vid_path1);
global shouldBreak;
shouldBreak = false;
global switch_vid;
switch_vid = true;

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A~=0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();      %neuron's center
gSiz = obj.options.gSiz;        % maximum size of a neuron

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% keep the log information
if ~save_img
    try
        log_file =  obj.P.log_file;
        flog = fopen(log_file, 'a');
        log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
        manual_intervention.before = obj.obj2struct();
        
        fprintf(flog, '[%s]\b', get_minute());
        fprintf(flog, 'Start manual interventions:\n');
    end
end

%% start viewing neurons
fig = figure('position', [100, 100, 1620, 720], 'color', [1 1 1]);
set(fig, 'Name', 'NeuraCleanse');
set(fig, 'KeyPressFcn', @key_press);
m=1;
s = size(obj.Cn);
h = s(1);
w = s(2);
su = sum(obj.A, 2);
nor = max(su);
su = reshape(full(su),h,w);


while and(m>=1, m<=length(ind))
    %% full-frame view
    ax1 = subplot(3,6,[1 2 7 8]); cla;
    sfps=reshape(full(obj.A(:, ind(m))),h,w);
%     plot_contours(obj.A(:,ind(m)), obj.Cn, 0.6, false, [], obj.Coor(ind(m)), 2);
    colormap(ax1, gray);
    Ed = edge(sfps, 'canny');
    B = imoverlay(obj.Cn,Ed,'r');
    imshow(B,'InitialMagnification',200);
    axis equal; axis off;
    
    tolerance = 0.05;
    if all(sum(sfps(1:floor(tolerance*w), :)) == 0) && all(sum(sfps(:, 1:floor(tolerance*h))) == 0) && all(sum(sfps(end-floor(tolerance*w):end, :)) == 0) && all(sum(sfps(:, end-floor(tolerance*h):end)) == 0)
        title('AlessandROI: This Neuron seems to be legit', 'Color', 'k');
    else
        title('AlessandROI: This Neuron is at the edge', 'Color', 'r');
    end
    axis equal; axis off;
    %% zoomed-in view
    ax2 = subplot(3,6,[5 6 11 12]); cla;
    dim = size(obj.A);
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    colormap(ax2, hot);
    %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    if ind_del(m)
        title(sprintf('Neuron %d/%d', ind(m), dim(2)), 'color', 'r');
    else
        title(sprintf('Neuron %d/%d', ind(m), dim(2)));
    end
    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end
    %% temporal components
    subplot(3,6,13:18);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m))), 'linewidth', 2); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r'); hold on;
        plot(t, obj.S(ind(m), :)*max(obj.A(:, ind(m))), 'magenta')
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    xlim([t(1), t(end)]);
    xlabel(str_xlabel);
    title('AlessandROI@https://mahechen.com:443 >> [K - Keep] [B - Back] [I - Inspect Video] [D - Delete] [R - Remove Rests] [E - End Session]', ind(m));
    %% save images
    videofig(vid.NumFrames, @(frm) playvideo(frm, vid));
    playvideo(1, vid);
    if save_img
        drawnow();
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
%     else
%         fprintf('Neuron %d, keep(k, default)/delete(d)/split(s)/trim(t)\n\t/trim cancel(tc)/delete all(r)/backward(b)/end(e):    ', ind(m));
    end
    while ~waitforbuttonpress
        pause(1)
    end
    if shouldBreak
        break;
    end         
end
    
if save_img
    cd(cur_cd);
else
    if ~isempty(ind(ind_trim))
        obj.A(:, ind(ind_trim)) = obj.A(:,ind(ind_trim)).*Amask(:, ind(ind_trim));
        try
            fprintf(flog, '\n\tFollowing neurons were trimmed:\n');
            ids_trimmed = ind(ind_trim);
            for m=1:length(ids_trimmed)
                fprintf(flog, '%2d, ', ids_trimmed(m));
            end
            fprintf(flog, '\n');
        end
    end
    
    if ~isempty(ind(ind_del))
        try
            fprintf(flog, '\tDeleting manually selected neurons:\n');
        end
        obj.delete(ind(ind_del));
    end
    %     obj.Coor = obj.get_contours(0.9);
    
    
    return;
end
try
    fprintf(flog, '[%s]\b', get_minute());
    fprintf(flog, 'Finished the manual intervention.\n');
    fprintf(flog, '[%s]\b', get_minute());
    if obj.options.save_intermediate
        manual_intervention.after = obj.obj2struct(); %#ok<STRNU>
        tmp_str = get_date();
        tmp_str=strrep(tmp_str, '-', '_');
        eval(sprintf('log_data.manual_%s = manual_intervention;', tmp_str));
        
        fprintf(flog, '\tThe results were saved as intermediate_results.manual%s\n\n', tmp_str);
    end
    fclose(flog);
end
    function playvideo(frame, vid)
        f = vid.read(frame);
        Ed = edge(sfps,'canny', 0.99);
        B = imoverlay(f,Ed,'r');
        imshow(B,'InitialMagnification',200);

    end
    function key_press(~, event)  %#ok, unused arguments
 		switch event.Key  %process shortcut keys
            case 'i'
                if switch_vid
                    vid = VideoReader(fullfile(vid_path, 'normvideo.avi'));
                    videofig(vid.NumFrames, @(frm) playvideo(frm, vid));
                    playvideo(1 ,vid)
                    switch_vid = false;
                else
                    vid = VideoReader(fullfile(vid_path, 'msvideo.avi'));
                    videofig(vid.NumFrames, @(frm) playvideo(frm, vid));
                    playvideo(1 ,vid)
                    switch_vid = true;
                end
            case 'd'
                ind_del(m) = true;
                m = m+1;
            case 'b'
                m = m-1;
            case 'r'
                ind_del(m:end) = true;
                shouldBreak = true;
            case 'k'
                ind_del(m) = false;
                m= m+1;
            case 's'
                try
                    subplot(336);
                    temp = imfreehand();
                    tmp_ind = temp.createMask();
                    tmpA = obj.A(:, ind(m));
                    obj.A(:, end+1) = tmpA.*tmp_ind(:);
                    obj.C(end+1, :) = obj.C(ind(m), :);
                    obj.A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                    obj.S(end+1, :) = obj.S(ind(m), :);
                    obj.C_raw(end+1, :) = obj.C_raw(ind(m), :);
                    obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(m), :);
                    k_ids = obj.P.k_ids;
                    obj.ids(end+1) = k_ids+1;   % assign an neuron id
                    obj.tags(end+1) = obj.tags(ind(m));
                    obj.P.k_ids = k_ids+1;
                    fprintf(flog, '\tSplit %d --> %d + %d\n', obj.ids(ind(m)),obj.ids(ind(m)), k_ids);
                catch
                    fprintf('the neuron was not split\n');
                end
            case 't'
                try
                    subplot(336);
                    temp = imfreehand();
                    tmp_ind = temp.createMask();
                    Amask(:, ind(m)) = tmp_ind(:);
                    ind_trim(m) = true;
                catch
                    fprintf('the neuron was not trimmed\n');
                end
            case 'c'
                Amask(:, ind(m)) = (obj.A(:, ind(m)) > 0);
                ind_trim(m) = false;
            case 'e'
                shouldBreak = true;
%             elseif ~isnan(str2double(temp))
%                 m = m + floor(str2double(temp));
%                 m = max(m, 1);
%                 m = min(m, length(ind));
%                 fprintf('jump to neuron %d / %d\n', m, length(ind));
%         else
%             m = m+1;
        end
    end
end