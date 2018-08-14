function SCycles=SizeTwoCycles(Cs)

% function SCycles=ShortCycles(Cs)

% Version: 3f
% Date: July 8, 2009

n=length(Cs);
SCycles=cell(0);

% Find directly revealed preferred items (x>y and y>x)
for x=1:n
    for y=(x+1):n
        if Cs(x,y)==1 && Cs(y,x)==1 % x>y and y>x
            SCycles{length(SCycles)+1}=[x y x]; % Store cycle
        end
    end
end

end