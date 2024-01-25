function A = AdjacencyMatrix(n)
A = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            A(i,j) = 1;
        end
    end
end
end

