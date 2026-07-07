function ref_points = gen_ref_points(ref_point_per_dimension, nobjs)
%GEN_REF_POINTS: generates the reference points on the unit hyperplane by
%enumerating all the disposition with fixed sum, i.e. all the point
%satisfying x_1 + x_2 + ... + x_nobjs = 1

ref_points = {};
gen_fix_sum_comb(ref_point_per_dimension,0,1,nobjs,zeros(1,nobjs));
ref_points = cell2mat(ref_points')./ref_point_per_dimension;

%% Nested Functions

function gen_fix_sum_comb(fixed_sum, curr_sum, depth, nobjs, curr_list)
    %GEN_FIX_SUM_COMB: recursively generates integer combinations having a 
    %   fixed sum, with fixed cardinality (specified by nobjs).
    %   depth should be set to 1 and curr_sum to 0 when calling this function
    %   curr_list should be set as zeros(1,nobjs)
    
    if depth == nobjs
        ret_num = fixed_sum - curr_sum;
        curr_list(depth) = ret_num;
        ref_points{end + 1} = curr_list;
        return
    else
        if fixed_sum > curr_sum
            for k = 0:(fixed_sum - curr_sum)
                curr_list(depth) = k;
                gen_fix_sum_comb(fixed_sum, curr_sum + k, depth + 1, nobjs, curr_list)
            end
        else
            curr_list(depth) = 0;
            gen_fix_sum_comb(fixed_sum, curr_sum, depth + 1, nobjs, curr_list)
    
        end
    end
end

end

