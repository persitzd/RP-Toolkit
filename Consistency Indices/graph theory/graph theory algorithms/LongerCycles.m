function LCycles=LongerCycles(Cs)

% function LCycles=LongerCycles(Cs)

% Version: 3f
% Date: July 8, 2009

n=length(Cs);
LCycles=cell(0);

% Find self-loops (x>x)
for zz=1:n
    if Cs(zz ,zz)>0
        LCycles{length(LCycles)+1}=[zz zz]; % Store cycle
    end
end

% Find directly revealed preferred items (x>y and y>x)
for x=1:n
    for y=(x+1):n
        if Cs(x,y)==1 && Cs(y,x)==1 % x>y and y>x
            LCycles{length(LCycles)+1}=[x y x]; % Store cycle
        end
    end
end

% Find 3 item cycles (x>y>z>x)
for x=1:n
    for y=(x+1):n
        for z=(y+1):n
            if Cs(x,y)==1 && Cs(y,z)==1 && Cs(z,x)==1
                LCycles{length(LCycles)+1}=[x y z x]; % Store cycle
            end
        end
    end
end

% Find 4 item cycles (x>y>z>a>x)
for x=1:n
    for y=(x+1):n
        for z=(y+1):n
            for a=(z+1):n
                if Cs(x,y)==1 && Cs(y,z)==1 && Cs(z,a)==1 && Cs(a,x)==1
                    stack=[x y z a];
                    [p,q,r]=dmperm(Cs(stack,stack)+eye(length(stack)));
                    teststack=max(diff(r));
                    if teststack<=1 || teststack==length(stack)
                        LCycles{length(LCycles)+1}=[x y z a x]; % Store cycle
                    end
                end
            end
        end
    end
end

% Find 5 item cycles (x>y>z>a>b>x)
for x=1:n
    for y=(x+1):n
        for z=(y+1):n
            for a=(z+1):n
                for b=(a+1):n
                    if Cs(x,y)==1 && Cs(y,z)==1 && Cs(z,a)==1 && Cs(a,b)==1 && Cs(b,x)==1
                        stack=[x y z a b];
                        [p,q,r]=dmperm(Cs(stack,stack)+eye(length(stack)));
                        teststack=max(diff(r));
                        if teststack<=1 || teststack==length(stack)
                            LCycles{length(LCycles)+1}=[x y z a b x]; % Store cycle
                        end
                    end
                end
            end
        end
    end
end

% Find 6 item cycles (x>y>z>a>b>c>x)
for x=1:n
    for y=(x+1):n
        for z=(y+1):n
            for a=(z+1):n
                for b=(a+1):n
                    for c=(b+1):n
                        if Cs(x,y)==1 && Cs(y,z)==1 && Cs(z,a)==1 && Cs(a,b)==1 && Cs(b,c)==1 && Cs(c,x)==1
                            stack=[x y z a b c];
                            [p,q,r]=dmperm(Cs(stack,stack)+eye(length(stack)));
                            teststack=max(diff(r));
                            if teststack<=1 || teststack==length(stack)
                               LCycles{length(LCycles)+1}=[x y z a b c x]; % Store cycle
                            end
                        end
                    end
                end
            end
        end
    end
end

end