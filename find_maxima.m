function [ new_max_fun, new_time ] = find_maxima( max_loc, i, N_points, s2, t2, Int_method, epsilon )
%FIND_MAXIMA Summary of this function goes here
%   Detailed explanation goes here
    index = max_loc(i);
    start_index = max(1, index - 2);
    end_index = min(N_points, index + 2);
    new_times = linspace(t2(start_index), t2(end_index), 100);
    s_sub = interp1(t2(start_index:end_index), s2(start_index:end_index), new_times, Int_method);
    [new_max_fun, max_loc_sub] = findpeaks(s_sub);
    new_time = new_times(max_loc_sub);
    [ new_max_fun, new_time, abs(s2(index)-new_max_fun)];
    if abs(s2(index)-new_max_fun) > epsilon
        [new_max_fun, new_time] = find_maxima(max_loc_sub, 1, 100, s_sub, new_times, Int_method, epsilon);
    end
end

