mice = ["038-01" "038-03" "038-04" "039-02" "039-04"];
sessions = ["day1" "day3" "day3am" "day5" "day5am" "imm" "pre" "recent" "training"];
analysis_file = "analysis_of_";
for mouse=1:5
    reg_data = load(fullfile("/media/jfadmin/LaCie/Mahe/baseline+day1/", mice(mouse), "Results/cellRegistered.mat"));
    dim = size(reg_data.cell_registered_struct.cell_to_index_map);
    analysis_struct = cell(flip(dim));
    for session=1:dim(2)
        data = load(fullfile("/media/jfadmin/LaCie/Mahe", sessions(session), int2str(mouse), "Result/data.mat"));
        S = data.st.S;
        trace_dim = size(S);
        for trace=1:dim(1)
            idx = reg_data.cell_registered_struct.cell_to_index_map(trace, session);
            if idx == 0
                analysis_struct{session, trace} = zeros(1, trace_dim(2));
            else
                analysis_struct{session, trace} = S(idx,:);
            end
        end
    end
    save(strcat("/media/jfadmin/LaCie/Mahe/baseline+day1/", analysis_file, mice(mouse), ".mat"), 'analysis_struct')
end