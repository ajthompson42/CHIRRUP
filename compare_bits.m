function propfound = compare_bits(input,output)

k1 = size(input,2);
k2 = size(output,2);

found = zeros(k1,1);
for i = 1:k1
    for j = 1:k2
        if any(output(:,j)~=input(:,i))
            continue
        end
        found(i,1) = j;
        break
    end
end

propfound = sum(found>0)/k1;