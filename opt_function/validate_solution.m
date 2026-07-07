function [flag,sat] = validate_solution(sat,Tmin_max,Tmin_max_cold,val_cases,T_op,T,T_cold,write_flag)

% Checks if the max and min temperature are inside the operative range for
% every component.
% The temperatures are expected to be in [°C] (Celsius).
% Returns 0 if the solution is valid, 1 if the solution is invalid.
% write_flag = 1 ---> results are saved inside "sat" struct

N_str_face = length(sat.geom.ext.face); % number of faces of the structure (equal to 6)
flag = 0; % 0 = valid, 1 = invalid

if isequal(val_cases,1)
    for l=1:sat.node.total_node
        % ### STRUCTURE NODES ###
        if l <= sat.node.ext.n_node % Hypothesis: structure 'ex' nodes always come first
            for k=1:length(sat.node.globe(l).face)
                face_id = sat.node.globe(l).face(k); 
                if (Tmin_max(l,1)<T_op(face_id,1) || Tmin_max(l,2)>T_op(face_id,2)) %check on hot case                        
                    flag=1;
                    break;
                end
            end
        else
            % ### ITEM NODES ###
            ID_T_op = sat.node.globe(l).ID_item - 1 + N_str_face; % id of T_op
            if (Tmin_max(l,1)<T_op(ID_T_op,1) || Tmin_max(l,2)>T_op(ID_T_op,2)) %check on hot case                        
                    flag=1;
                    break;
            end            
        end     

        if isequal(flag,1)
            break;
        elseif write_flag == 1
            sat.node.globe(l).Temperature_t=T(:,l);
        end
    end
end
if isequal(val_cases,2)
    for l=1:sat.node.total_node
        % ### STRUCTURE NODES ###
        if l <= sat.node.ext.n_node % Hypothesis: structure 'ex' nodes always come first
            for k=1:length(sat.node.globe(l).face)
                face_id = sat.node.globe(l).face(k); 
                if (Tmin_max(l,1)<T_op(face_id,1) ||...
                        Tmin_max(l,2)>T_op(face_id,2))||... %check on hot case
                        (Tmin_max_cold(l,1)<T_op(face_id,1) ||...
                        Tmin_max_cold(l,2)>T_op(face_id,2)) %check on cold case                       
                    flag=1;
                    break;
                end
            end
        else
            % ### ITEM NODES ###
            ID_T_op = sat.node.globe(l).ID_item - 1 + N_str_face; % id of T_op
            if (Tmin_max(l,1)<T_op(ID_T_op,1) ||...
                    Tmin_max(l,2)>T_op(ID_T_op,2))||... %check on hot case
                    (Tmin_max_cold(l,1)<T_op(ID_T_op,1) ||...
                    Tmin_max_cold(l,2)>T_op(ID_T_op,2)) %check on cold case                       
                flag=1;
                break;
            end            
        end

        if isequal(flag,1)
            break;
        elseif write_flag == 1
            sat.node.globe(l).Temperature_t=T(:,l);
            sat.node.globe(l).Temperature_t_cold=T_cold(:,l);
        end
    end
end

end

