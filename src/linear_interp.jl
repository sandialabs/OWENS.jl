function [ynew] = linear_interp(x_array,y_array,xnew)
ynew = zeros(size(xnew));
for i = 1:length(xnew)
    ynew(i) = sub_linear_interp(x_array,y_array,xnew(i));
end
end

function [ynew] = sub_linear_interp(x_array,y_array,xnew)


i = 1;
min_x = min(x_array);
max_x = max(x_array);
if xnew<min_x # extrapolate linearly at max and min
#     ynew = y_array(min_x_idx);
    i = 1;
elseif xnew>max_x
    #     ynew = y_array(max_x_idx);
    i = length(x_array)-1;
else
    for j = 1:length(x_array)-1
        if x_array(j)<=xnew && x_array(j+1)>=xnew
            i = j;
            break            
        end
    end    
end

fraction = (xnew - x_array(i))/(x_array(i+1) - x_array(i));
ynew = y_array(i) + fraction * (y_array(i+1) - y_array(i));

end
