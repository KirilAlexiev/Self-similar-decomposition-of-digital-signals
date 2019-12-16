function [indicator, level]=find_similar(a_in, indicator, level, extr_arrS, t2, s2, Int_method)
     %============ Calculate similarity ===================================
     len_a_in = length(a_in);
     s = zeros(100,len_a_in-1);
     corr_data = zeros(1,len_a_in);
     for i = 1:len_a_in
         ind1 = extr_arrS(a_in(i),1);
         ind2 = extr_arrS(a_in(i)+1,1);
         dt = (t2(ind2)-t2(ind1))/99;
         s(:,i) = interp1(t2(ind1:ind2), s2(ind1:ind2), t2(ind1):dt:t2(ind2), Int_method);
         corr_data(i)=corr(s(:,1),s(:,i));
     end
     if level == 0
        indicator(a_in(1)) = 1;
     end
     equal_data_ind = find(corr_data > 0.98);
     len_eq = length(equal_data_ind);
     if len_eq > 1
        equal_data = a_in(equal_data_ind);
        level = level +1;
        indicator(equal_data) = level;
        a_in_next = equal_data + 1;   %search for the next  el+one;
        % check for existance of elements
        not_ind = [];
        for i = 1:length(a_in_next)
             nonzero_el = find(a_in==a_in_next(i));
             if nonzero_el
                not_ind = [not_ind, i];
             end
        end
        a_in_next = a_in_next(not_ind);
        ll = length(a_in_next);
        switch ll
            case 0 % nothing to do
            case 1
                 indicator(a_in(1)) = 1;
            otherwise
                 for i = 1: length(a_in_next)-1
                    [indicator, level]=find_similar(a_in_next, indicator, level, extr_arrS, t2, s2, Int_method);
                 end
        end
        level = level - 1;
     end
end