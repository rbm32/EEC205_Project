function tf = is_double_single_4d_matrix(x)
    tf =  (...
            strcmp(string(class(x)),"double") ...
            || strcmp(string(class(x)),"single")...
           )...
          && (numel(size(x)) == 4);
end