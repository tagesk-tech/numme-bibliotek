%simplex

% 

f  = [1; 2; 3];          % kostnadsvektor för minimering
A  = [-2 -1 -1;          % ≥-villkor skrivs om som ≤ med minus
      -1 -3 -2;
      -1 -1 -4];
b  = [-5; -7; -6];

[x, index, y] = simplex(A, b, f); 

function [x, index, y, val] = simplex(B, b, f)

n = numel(b);
I = eye(n);
A = [B, I]; % builds the constriant with slack variables
c = [f; zeros(n, 1)];
[n, m] = size(A);
index = [n+1:m]; % gets the starting index values
non_index = [1:n]; % stores the locked variables

for i = 1:1 % starts the indexing loop
    A_b = A(:, index); % gets the first solution
    A_v = A(:, non_index)
    c_v = c(non_index);
    c_b = c(index); 
    y = A_b'\c_b; % calculates the new y
    y
    y'* A_v
    r_v = c_v' - y' * A_v;

    [min_rv, idx_add] = min(r_v);
    if min_rv >= 0
        break;
    end

    a_k = A_v(:, idx_add);
    a_bar = A_b\a_k; %calculates the a_bar
    b_bar = A_b\b; 
    t = b_bar./a_bar;
    [val, idx_remove] = mint(t); %tittar på dem indexen som faktiskt är med
    %now we are ready to switch in this case
    idx1 = index(idx_remove);
    idx2 = non_index(idx_add);
    index;
    non_index;
    
    index(idx_remove) = idx2;
    non_index(idx_add) = idx1;

end

m = size(A,2);
x = zeros(m,1);
x(index) = A_b \ b;    % baslösningen

y = A_b' \ c(index);



% when we have the simplex we can solve the dual in this case


end 




function [val, idx] = mint(t)
    val = Inf; idx = NaN;
    for i = 1:numel(t)
        if t(i) >= 0 && t(i) < val
            val = t(i);
            idx = i;
        end
    end
end