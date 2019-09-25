addpath("utils")
r = 0
l = 0
re = 1

EbN0=10
trials = 1

  
prop=1  
mvalues=[6]
pvalues=[5]

for a = 1:size(mvalues, 2)
    m=mvalues(a);
    
    for b=1:size(pvalues,2)
        p=pvalues(b);
        i=1;
        disp('m')
        disp(m)
        disp('p')
        disp(p)
        output=[];
        K=[] ;
        time=[];
        prop=1;

        
        while prop > 0.05
            K=[K, i];
            disp(i)
            [ave_time, prop] = run(r, l, re, m, p, EbN0, i, trials);
            output = [output, prop];
            time = [time, ave_time];
            disp(prop)
            i=i+1;
        end
        
        filename=strcat("tests/B", num2str(B),'r', num2str(r),'l', num2str(l),'r', num2str(r),'m', num2str(m),'p', num2str(p), 'trials', num2str(trials))
        save(filename, "K", "output", "time");
    end
end



x=strcat("tests/B", num2str(B),'r', num2str(r),'l', num2str(l),'r', num2str(r),'m', num2str(m),'p', num2str(p), 'trials', num2str(trials))
x(1)
plot(K, output)
xlabel('K')
ylabel('prop found')