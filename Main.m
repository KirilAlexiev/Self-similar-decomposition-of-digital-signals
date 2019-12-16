close all
% Signal parameters
%f = 31; %20   f = (nn-1)/(w/(2*pi))     % sampling rate
nn = 69; %35; %60                       % number of sampling points
n_max = floor((nn-1)/2);
discrets_ = 512;                         % number of discretization levels
t = (0:(nn-1))/(nn-1);                   % Input vector of sampling moments
Time_interval = 1;                       % signal longevity
w = 9*pi;                                % 2*pi Observed signal

epsilon = 0.00001;                       % precision

%============= Examples - Choose one and comment others ===============
%============= Example 1 
%
s_anal = sin(w*t);

%============= Example 2 
%s_anal = t.*sin(w*t);                       

%============= Example 3 
%s_anal = sin(w*(t.*t));                      

%============= Example 4 
%s_anal = sin(w*t) .* (1 + w*t / 2);

%============= Example 5 
%s_anal = sin(w*t) + sin(w*t/2 + 0.75);      

%============= Example 6 
%s_anal = sin(w*t) + w*t * 0.2;              

%============= Example 7 (counterexample)
%s_anal = sin(w*t) - (w*t.^2);               

%============= Example 8 
%s_anal = diric(w*t, 5);

%============= Example 9 
%s_anal = sawtooth(w*t-0.01);                

%============= Example 10 
%w = 9.01*pi;
%s_anal = square(w*t-0.01);                   

%============= Example 11 
%s_anal = sinc(w*t-pi);

%============= Common part ========================
N_points = 500;

%============= Example 12 Special input for ECG signal =================
% load('ECG','dataOut');
% s_anal = dataOut(10:end-2000);
% N_points = length(s_anal)
% t = 0:1/360:(N_points-1)/360;   % pulse beat in minute
% N_points = 5000;
% Time_interval = t(end);
%========================== end input data for examples ===================
 

%============ Filtration if necessary ==================
%???
%============ Initialization - Choose interpolation method ================
%Int_method = 'spline';
Int_method = 'cubic';
%Int_method = 'linear';

s_max = max(s_anal);
s_min = min(s_anal);
max_spread = s_max-s_min;
s_discr = floor(((s_anal-s_min)./max_spread)*discrets_); % Discrete signal
figure, plot(t,s_discr,'b-x');
t2 = ((0:(N_points-1))/(N_points-1)).*Time_interval;      % New sampling points
s2 = interp1(t,s_discr,t2,Int_method);                    % New signal points
hold on, plot(t2,s2,'-o');

%============ Find Extrema =============================
s_max = max(s2);
s_min = min(s2);
max_spread = s_max-s_min;
s_discr = floor(((s2-s_min)./max_spread)*discrets_); % Discrete signal
s2 = s_discr;

[max_fun,max_loc] = findpeaks(s2);
for i = 1:length(max_loc)
    [new_value, new_time] = find_maxima(max_loc, i, N_points, s2, t2, Int_method, epsilon);
    t2 = sort([t2, new_time]);
    index = max_loc(i);
    if new_time > t2(index)
        s2 = [s2(1:index), new_value, s2(index+1:end)];
    else
        s2 = [s2(1:index-1), new_value, s2(index:end)];
    end
    max_loc = max_loc+1;
end
[min_fun,min_loc] = findpeaks(-s2);
for i = 1:length(min_loc)
    [new_value, new_time] = find_maxima(min_loc, i, N_points, -s2, t2, Int_method, epsilon);
    new_value = -new_value;
    t2 = sort([t2, new_time]);
    index = min_loc(i);
    if new_time > t2(index)
        s2 = [s2(1:index), new_value, s2(index+1:end)];
    else
        s2 = [s2(1:index-1), new_value, s2(index:end)];
    end
    min_loc = min_loc+1;
end
[max_fun,max_loc] = findpeaks(s2);
[min_fun,min_loc] = findpeaks(-s2);
min_fun = -min_fun;
extr_arr = [max_loc', max_fun'; min_loc', min_fun'];
extr_arrS = sortrows(extr_arr,1);
[sN, sM] = size(extr_arrS);              % sN - Number of extrema
N_segments = sN - 1;                     % N_segment - Number of segments under consideration
%============= Main procedure ===========================
indicator = zeros(1, N_segments);
All_shapes = zeros(N_segments,2);
Originality = zeros(1,N_segments);
for i = 1: N_segments
    a_in = i:N_segments; level = 0;
    nonzero_el = find(indicator(a_in)>0);
    if nonzero_el 
       a_in(nonzero_el) = [];
    end 
    ll = length(a_in);
    switch ll
        case 0 % nothing to do
        case 1
           indicator(a_in(1)) = 1;
        otherwise
           [indicator, level]=find_similar(a_in, indicator, level, extr_arrS, t2, s2, Int_method);
    end
end
%=================== print segments ======================
figure,
col_pal = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0];
i_col = 1;
plot(t2(1:(extr_arrS(1,1)-1)),s2(1:(extr_arrS(1,1)-1)),'o','Color',col_pal(i_col,:)); hold on;
for i = 1:N_segments
    if indicator(i) == 1 
       i_col = i_col + 1;
       if i_col == 8 
          i_col = 1; 
       end
       plot(t2(extr_arrS(i,1):extr_arrS(i+1,1)-1),...
            s2(extr_arrS(i,1):extr_arrS(i+1,1)-1),'o','Color',col_pal(i_col,:)); hold on;
    end;
    plot(t2(extr_arrS(i,1):extr_arrS(i+1,1)-1),...
         s2(extr_arrS(i,1):extr_arrS(i+1,1)-1),'o','Color',col_pal(i_col,:)); hold on;
end
i_col = i_col + 1;
plot(t2(extr_arrS(sN,1):end),...
         s2(extr_arrS(sN,1):end),'o','Color',col_pal(i_col,:));
ax_sc = axis;     
grid on;
ax = gca;
ax.FontSize = 16;
% ================== print primitives ===================== 
i_counter = 0;
for i = 1:N_segments
    if (indicator(i)==1)
        i_counter = i_counter + 1;
        All_shapes(i_counter,1) = extr_arrS(i,1);
    end
    All_shapes(i_counter,2) = extr_arrS(i+1,1)-1;
end
Originality(1:i_counter) = 1;
for i = 1:i_counter-1
    for j = i+1:i_counter
        ind11 = All_shapes(i,1);
        ind12 = All_shapes(i,2);
        dt1 = (t2(ind12)-t2(ind11))/99;
        ind21 = All_shapes(j,1);
        ind22 = All_shapes(j,2);
        dt2 = (t2(ind22)-t2(ind21))/99;
        ss1 = interp1(t2(ind11:ind12), s2(ind11:ind12), t2(ind11):dt1:t2(ind12), Int_method);
        ss2 = interp1(t2(ind21:ind22), s2(ind21:ind22), t2(ind21):dt2:t2(ind22), Int_method);
        corr_=corr(ss1',ss2');
        if corr_>0.96
           Originality(j) = 0;
        end
    end
end

figure, 
col_pal = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0];
i_col = 1;
plot(t2(1:(extr_arrS(1,1)-1)),s2(1:(extr_arrS(1,1)-1)),'o','Color',col_pal(i_col,:)); hold on;
for i = 1:i_counter
    if Originality(i) == 1
       i_col = i_col + 1;
       if i_col == 8 
          i_col = 1; 
       end
       plot(t2(All_shapes(i,1):All_shapes(i,2)),...
            s2(All_shapes(i,1):All_shapes(i,2)),'o','Color',col_pal(i_col,:)); hold on;
    end;
end
i_col = i_col + 1;
if i_col == 8 
          i_col = 1; 
end
plot(t2(extr_arrS(sN,1):end),...
         s2(extr_arrS(sN,1):end),'o','Color',col_pal(i_col,:));
axis(ax_sc);     
grid on
ax = gca;
ax.FontSize = 16;
        
       
