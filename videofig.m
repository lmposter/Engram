function [fig_handle, axes_handle, scroll_bar_handles, scroll_func] = ...
	videofig(num_frames, redraw_func, play_fps, big_scroll, ...
	key_func, varargin)
	
	%default parameter values
	if nargin < 3 || isempty(play_fps), play_fps = 30; end  %play speed (frames per second)
	if nargin < 4 || isempty(big_scroll), big_scroll = 30; end  %page-up and page-down advance, in frames
	if nargin < 5, key_func = []; end
	
	%check arguments
	check_int_scalar(num_frames);
	check_callback(redraw_func);
	check_int_scalar(play_fps);
	check_int_scalar(big_scroll);
	check_callback(key_func);

	click = 0;
	f = 1;  %current frame
    
    %initialize figure
	fig_handle = gcf;
    set(fig_handle, 'MenuBar','none', 'Units','norm', ...
		'WindowButtonDownFcn',@button_down, 'WindowButtonUpFcn',@button_up, ...
		'WindowButtonMotionFcn', @on_click,  ...
		'Interruptible','off', 'BusyAction','cancel', varargin{:});
	
	
	%axes for scroll bar
	scroll_axes_handle = axes('Parent',fig_handle, 'Position',[0.12 0.11 0.78 0.03], ...
		'Visible','off', 'Units', 'normalized');
	axis([0 1 0 1]);
	axis off
	
	%scroll bar
	scroll_bar_width = max(1 / num_frames, 0.01);
	scroll_handle = patch([0 1 1 0] * scroll_bar_width, [0 0 1 1], [.8 .8 .8], ...
		'Parent',scroll_axes_handle, 'EdgeColor','b', 'ButtonDownFcn', @on_click);
	
	%timer to play video
	play_timer = timer('TimerFcn',@play_timer_callback, 'ExecutionMode','fixedRate');
	
	%main drawing axes for video display
	axes_handle = axes('Position',[0.2075 0.38 0.62 0.62]);
	
	%return handles
	scroll_bar_handles = [scroll_axes_handle; scroll_handle];
	scroll_func = @scroll;
	
	
	
% 	function key_press(src, event)  %#ok, unused arguments
% 		switch event.Key,  %process shortcut keys
% 		case 'leftarrow',
% 			scroll(f - 1);
% 		case 'rightarrow',
% 			scroll(f + 1);
% 		case 'pageup',
% 			if f - big_scroll < 1,  %scrolling before frame 1, stop at frame 1
% 				scroll(1);
% 			else
% 				scroll(f - big_scroll);
% 			end
% 		case 'pagedown',
% 			if f + big_scroll > num_frames,  %scrolling after last frame
% 				scroll(num_frames);
% 			else
% 				scroll(f + big_scroll);
% 			end
% 		case 'home',
% 			scroll(1);
% 		case 'end',
% 			scroll(num_frames);
% 		case 'return',
% 			play(1/play_fps)
% 		case 'backspace',
% 			play(5/play_fps)
% 		otherwise,
% 			if ~isempty(key_func),
% 				key_func(event.Key);  %#ok, call custom key handler
% 			end
% 		end
% 	end
	
	%mouse handler
	function button_down(src, event)  %#ok, unused arguments
		set(src,'Units','norm')
		click_pos = get(src, 'CurrentPoint');
        if (0.12 <= click_pos(1)) && (click_pos(1) <= 0.88)
            if (0.1 <= click_pos(2)) && (click_pos(2) <= 0.13)  %only trigger if the scrollbar was clicked
                click = 1;
                on_click([],[]);
            end
        end
	end

	function button_up(src, event)  %#ok, unused arguments
		click = 0;
	end

	function on_click(src, event)  %#ok, unused arguments
		if click == 0, return; end
		
		%get x-coordinate of click
		set(fig_handle, 'Units', 'normalized');
		click_point = get(fig_handle, 'CurrentPoint');
		set(fig_handle, 'Units', 'pixels');
		x = click_point(1);
        if x <= 0.51
            x = x - 0.12 * (1 - ((x - 0.12)/0.51));
        else
            x = x + 0.12 * ((x - 0.51)/0.51);
        end

		%get corresponding frame number
		new_f = floor(1 + x * num_frames);
		
		if new_f < 1 || new_f > num_frames, return; end  %outside valid range
		
		if new_f ~= f,  %don't redraw if the frame is the same (to prevent delays)
			scroll(new_f);
		end
	end

	function play(period)
		%toggle between stoping and starting the "play video" timer
		if strcmp(get(play_timer,'Running'), 'off'),
			set(play_timer, 'Period', period);
			start(play_timer);
		else
			stop(play_timer);
		end
	end
	function play_timer_callback(src, event)  %#ok
		%executed at each timer period, when playing the video
		if f < num_frames,
			scroll(f + 1);
		elseif strcmp(get(play_timer,'Running'), 'on'),
			stop(play_timer);  %stop the timer if the end is reached
		end
	end

	function scroll(new_f)
		if nargin == 1,  %scroll to another position (new_f)
			if new_f < 1 || new_f > num_frames,
				return
			end
			f = new_f;
		end
		
		%convert frame number to appropriate x-coordinate of scroll bar
		scroll_x = (f - 1) / num_frames;
		
		%move scroll bar to new position
		set(scroll_handle, 'XData', scroll_x + [0 1 1 0] * scroll_bar_width);
		
		%set to the right axes and call the custom redraw function
		set(fig_handle, 'CurrentAxes', axes_handle);
		redraw_func(f);
		
		%used to be "drawnow", but when called rapidly and the CPU is busy
		%it didn't let Matlab process events properly (ie, close figure).
		pause(0.001)
	end
	
	%convenience functions for argument checks
	function check_int_scalar(a)
		assert(isnumeric(a) && isscalar(a) && isfinite(a) && a == round(a), ...
			[upper(inputname(1)) ' must be a scalar integer number.']);
	end
	function check_callback(a)
		assert(isempty(a) || strcmp(class(a), 'function_handle'), ...
			[upper(inputname(1)) ' must be a valid function handle.'])
	end
end