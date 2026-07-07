function [U,U_T_low_hot,U_T_high_hot,U_T_low_cold,U_T_high_cold,U_H_hot,U_H_cold,U_L,...
    U_T_low_hot_f,U_T_high_hot_f,U_T_low_cold_f,U_T_high_cold_f,U_H_hot_f,U_H_cold_f,U_L_f]=...
    create_t_utility(sat,opt_st,Tmin_max_st,heat_sum_h,heat_sum_c,Vol_link_sum,val_cases,ROB_STRUCT)
% compute utility of solutions, considering every objective

N_str_face = length(sat.geom.ext.face); % number of faces of the structure (equal to 6)
T_op=opt_st.T_op_item;
T_op_all=zeros(size(T_op,1),4); % 1 = high-hot/default; 2 = low-hot/default; 3 = high-cold; 4 = low-cold;
gap_diff_hot=zeros(size(Tmin_max_st,2),2); % hot/default
gap_diff_cold=zeros(size(Tmin_max_st,2),2); % cold

%% Temperature gap computation

for i=1:length(Tmin_max_st) % for every solution...

    % T_op_all initialization
    T_op_all(:,1)=180; % +Inf initialization
    T_op_all(:,2)=-180; % -Inf initialization
    if isequal(val_cases,2)
        T_op_all(:,3)=180; % +Inf initialization
        T_op_all(:,4)=-180; % -Inf initialization
    end

    for j=1:sat.node.total_node  % for every node...
        
        % ### STRUCTURE NODES ###
        if j <= sat.node.ext.n_node % Hypothesis: structure 'ex' nodes always come first
            for k=1:length(sat.node.globe(j).face)
                face_id = sat.node.globe(j).face(k);
                if Tmin_max_st(i).data(j,1)<T_op_all(face_id,1) % high temperature side, in the hot/default case 
                    T_op_all(face_id,1)=Tmin_max_st(i).data(j,1);
                end
                if Tmin_max_st(i).data(j,2)>T_op_all(face_id,2) % low temperature side, in the hot/default case 
                    T_op_all(face_id,2)=Tmin_max_st(i).data(j,2);
                end
                if isequal(val_cases,2)
                    if Tmin_max_st(i).data_cold(j,1)<T_op_all(face_id,3) % high temperature side, in the cold case 
                        T_op_all(face_id,3)=Tmin_max_st(i).data_cold(j,1);
                    end
                    if Tmin_max_st(i).data_cold(j,2)>T_op_all(face_id,4) % low temperature side, in the cold case 
                        T_op_all(face_id,4)=Tmin_max_st(i).data_cold(j,2);
                    end
                end

            end
        else
            % ### ITEM NODES ###
            ID_T_op = sat.node.globe(j).ID_item - 1 + N_str_face; % id of T_op_all
            if Tmin_max_st(i).data(j,1)<T_op_all(ID_T_op,1) % high temperature side, in the hot/default case 
                T_op_all(ID_T_op,1)=Tmin_max_st(i).data(j,1);
            end
            if Tmin_max_st(i).data(j,2)>T_op_all(ID_T_op,2) % low temperature side, in the hot/default case 
                T_op_all(ID_T_op,2)=Tmin_max_st(i).data(j,2);
            end
            if isequal(val_cases,2)
                if Tmin_max_st(i).data_cold(j,1)<T_op_all(ID_T_op,3) % high temperature side, in the cold case 
                    T_op_all(ID_T_op,3)=Tmin_max_st(i).data_cold(j,1);
                end
                if Tmin_max_st(i).data_cold(j,2)>T_op_all(ID_T_op,4) % low temperature side, in the cold case 
                    T_op_all(ID_T_op,4)=Tmin_max_st(i).data_cold(j,2);
                end
            end
        end       

    end

    gap_hot=abs(T_op_all(:,1:2)-T_op); % temperature difference in the hot/default case (1 = T_op_min; 2 = T_op_max)
    gap_diff_hot(i,:)=[min(gap_hot(:,1)),min(gap_hot(:,2))]; % minimum gap in the low (1) and high (2) side, in the hot/default case 
    if isequal(val_cases,2)
        gap_cold=abs(T_op_all(:,3:4)-T_op); % temperature difference in the cold case (1 = T_op_min; 2 = T_op_max)
        gap_diff_cold(i,:)=[min(gap_cold(:,1)),min(gap_cold(:,2))]; % minimum gap in the low (1) and high (2) side, in the cold case 
    end
end

%% Utility computation

opt_st.obj_weight = [opt_st.obj_weight(1), opt_st.obj_weight]; % adding one utility objective (both low and high side)
opt_st.objective_active = [opt_st.objective_active(1), opt_st.objective_active];
if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)
    U_T_low_hot_f=ROB_STRUCT.U_T_low_hot_f;
    U_T_high_hot_f=ROB_STRUCT.U_T_high_hot_f;
else
    U_T_low_hot_f=fit([0 max(gap_diff_hot(:,1))]',[0 1]','poly1'); % scale function: from 0 to minimum gap in the low side (hot/default case)
    U_T_high_hot_f=fit([0 max(gap_diff_hot(:,2))]',[0 1]','poly1'); % scale function: from 0 to minimum gap in the high side (hot/default case)
end
U_T_low_hot = U_T_low_hot_f(gap_diff_hot(:,1));
U_T_high_hot = U_T_high_hot_f(gap_diff_hot(:,2));
utility_data=[U_T_low_hot U_T_high_hot];

U_T_low_cold = [];
U_T_low_cold_f = [];
U_T_high_cold = [];
U_T_high_cold_f = [];
if isequal(val_cases,2)
    opt_st.obj_weight = [opt_st.obj_weight(1:2), opt_st.obj_weight(3), opt_st.obj_weight(3:end)]; % adding one utility objective (both low and high side)
    opt_st.objective_active = [opt_st.objective_active(1:2), opt_st.objective_active(3), opt_st.objective_active(3:end)];
    if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)
        U_T_low_cold_f=ROB_STRUCT.U_T_low_cold_f;
        U_T_high_cold_f=ROB_STRUCT.U_T_high_cold_f;
    else
        U_T_low_cold_f=fit([0 max(gap_diff_cold(:,1))]',[0 1]','poly1'); % scale function: from 0 to minimum gap in the low side (hot/default case)
        U_T_high_cold_f=fit([0 max(gap_diff_cold(:,2))]',[0 1]','poly1'); % scale function: from 0 to minimum gap in the high side (hot/default case)
    end
    U_T_low_cold = U_T_low_cold_f(gap_diff_cold(:,1));
    U_T_high_cold = U_T_high_cold_f(gap_diff_cold(:,2));
    utility_data=[utility_data U_T_low_cold U_T_high_cold];
end

U_H_hot = [];
U_H_hot_f = [];
if opt_st.N_heater_hot > 0
    if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)
        U_H_hot_f=ROB_STRUCT.U_H_hot_f; 
    else
        U_H_hot_f=fit([min(heat_sum_h) max(heat_sum_h)]',[1 0]','poly1'); % scale function: from minimum total heater power to maximum total heater power (hot/default case)
        if min(heat_sum_h) == max(heat_sum_h) && min(heat_sum_h) == 0 % degenerate case
            U_H_hot_f = @(x) 1*ones(size(x));
        end
    end
    U_H_hot = U_H_hot_f(heat_sum_h);
    utility_data=[utility_data U_H_hot];
end

U_H_cold = [];
U_H_cold_f = [];
if opt_st.N_heater_cold > 0
    if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)
        U_H_cold_f=ROB_STRUCT.U_H_cold_f;
    else
        U_H_cold_f=fit([min(heat_sum_c) max(heat_sum_c)]',[1 0]','poly1'); % scale function: from minimum total heater power to maximum total heater power (cold case)
        if min(heat_sum_c) == max(heat_sum_c) && min(heat_sum_c) == 0 % degenerate case
            U_H_cold_f = @(x) 1*ones(size(x));
        end
    end    
    U_H_cold = U_H_cold_f(heat_sum_c);
    utility_data=[utility_data U_H_cold];
end

U_L = [];
U_L_f = [];
if opt_st.N_link > 0
    if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)
        U_L_f=ROB_STRUCT.U_L_f;
    else
         U_L_f=fit([min(Vol_link_sum) max(Vol_link_sum)]',[1 0]','poly1'); % scale function: from minimum total heater power to maximum total heater power (cold case)
         if min(Vol_link_sum) == max(Vol_link_sum) && min(Vol_link_sum) == 0 % degenerate case
             U_L_f = @(x) 1*ones(size(x));
         end
    end
    U_L = U_L_f(Vol_link_sum);
    utility_data=[utility_data U_L];
end

active_obj_weight = repmat(opt_st.obj_weight(opt_st.objective_active == 1), size(utility_data,1), 1);
% weighted average of utility_data using opt_st.obj_weight:
% utility_data is n_solution x n_utility_objectives
U = sum(utility_data.*(active_obj_weight),2)./(sum(opt_st.obj_weight(opt_st.objective_active == 1))); 

if any(isnan(U))
    disp("U NaN")
end

end