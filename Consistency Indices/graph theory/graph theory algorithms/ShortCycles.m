function SCycles=ShortCycles(Cs)

% function SCycles=ShortCycles(Cs)

% Version: 3f
% Date: July 8, 2009

n=length(Cs);
SCycles=cell(0);

% Find self-loops (x>x)
for zz=1:n
    if Cs(zz ,zz)>0
        SCycles{length(SCycles)+1}=[zz zz]; % Store cycle
    end
end

% Find directly revealed preferred items (x>y and y>x)
for x=1:n
    for y=(x+1):n
        if Cs(x,y)==1 && Cs(y,x)==1 % x>y and y>x
            SCycles{length(SCycles)+1}=[x y x]; % Store cycle
        end
    end
end

% Find 3 item cycles (x>y>z>x)
for x=1:n
    for y=(x+1):n
        for z=(y+1):n
            if Cs(x,y)==1 && Cs(y,z)==1 && Cs(z,x)==1
                SCycles{length(SCycles)+1}=[x y z x]; % Store cycle
            end
        end
    end
end

end