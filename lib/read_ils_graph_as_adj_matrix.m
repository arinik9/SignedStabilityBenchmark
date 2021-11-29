% https://fr.mathworks.com/matlabcentral/answers/355383-how-to-construct-an-undirected-graph-by-reading-edge-by-edge-from-a-text-file
function A = read_ils_graph_as_adj_matrix ( filepath )
    fid = fopen(filepath);
    C = textscan(fid,'%s%s%s');
    fclose(fid);

    m = size(C{1},1)-1;
    n = str2double(C{1}(1));

    A = zeros(n,n);

    V1 = str2double(C{1})+1;
    V2 = str2double(C{2})+1;
    W = str2double(C{3});

    for i = 2:(m+1)
        A(V1(i),V2(i)) = W(i);
        A(V2(i),V1(i)) = W(i);
    end
end
