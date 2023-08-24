mice = ["038-01" "038-03" "038-04" "039-02" "039-04"];
ref = "recent";
sess = ["day3" "day5" "imm" "pre" "training"];
analysis_file = "analysis_of_";
for mouse=5:5
    reg_data = load(fullfile("/media/jfadmin/LaCie/Mahe/ensemble_cellreg/", mice(mouse), "Results/cellRegistered.mat"));
    dim = size(reg_data.cell_registered_struct.cell_to_index_map);
    analysis_struct = cell(2, dim(1));
    for i=1:5
        sessions = [ref, sess(i)];
        for session=1:2
            data = load(fullfile("/media/jfadmin/LaCie/Mahe", sessions(session), int2str(mouse), "Result/data.mat"));
            S = data.st.C;
    %         S(S < 0) = 0;
            trace_dim = size(S);
            if session == 1
                id = 5;
            else
                if i == 5
                    id = 6;
                else
                    id = i;
                end
            end
            for trace=1:dim(1)
                
                idx = reg_data.cell_registered_struct.cell_to_index_map(trace, id);
                if idx == 0
                    analysis_struct{session, trace} = [];
                    disp(['NaN'])
                else
                    analysis_struct{session, trace} = S(idx,:);
                end
            end
        end
        emptyColumns = any(cellfun(@isempty, analysis_struct));
        analysis_struct = analysis_struct(:,~emptyColumns);
        save(strcat("/media/jfadmin/LaCie/Mahe/ensemble/", analysis_file, mice(mouse), sess(i), ".mat"), 'analysis_struct')
    end 
end