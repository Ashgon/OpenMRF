function fun_path = pulseq_get_path(fun_name)
    temp = strsplit(which(fun_name), filesep);
    for j=1:numel(temp)-1
        if j==1
            fun_path = [ temp{j} '/'];
        else
            fun_path = [ fun_path temp{j} '/'];
        end
    end
end