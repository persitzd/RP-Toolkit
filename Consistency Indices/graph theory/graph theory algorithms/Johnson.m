function [JCycles, R_Components]=Johnson(Cs)

% function JCycles=Johnson(Cs)
% List all cycles formed by unique (distinct) preference pairs.
% Comments prefixed with 'J:' are pseudocode from Johnson (1975).
% Modifications to Johnson's algorithm are noted.  Primarily, only
% cycles without subcycles are recorded.

% Version: 3f
% Date: July 8, 2009

global AK B blocked s items JCycles Ds Fs n

% Initialize global variables
Fs=Cs;
Ds=Cs;
n=length(Fs);
B=zeros(n,n);
blocked=zeros(1,n);
JCycles=cell(0);

% % Find self-loops (x>x)
% for zz=1:n
%     if Cs(zz,zz)==1
%         JCycles{length(JCycles)+1}=[zz zz]; % Store cycle
%         Cs(zz,zz)=0;
%     end
% end
% 
% % Find directly revealed preferred items (x>y and y>x)
% for x=1:n
%     for y=(x+1):n
%         if Cs(x,y)==1 && Cs(y,x)==1 % x>y and y>x
%             JCycles{length(JCycles)+1}=[x y x]; % Store cycle
%             Cs(x,y)=0;
%             Cs(y,x)=0;
%         end
%     end
% end

    [pp,qq,rr]=dmperm(Ds+eye(n));
    diffr=diff(rr); 
    for aa=1:length(diffr)
        R_Components(aa,:)=[sort(pp(rr(aa):rr(aa+1)-1)) zeros(1,n-diffr(aa))];
    end
    R_Components=sortrows(R_Components,1); % Sort by least vertex
    R_Components=R_Components(R_Components(:,2)>0,:); % Discard singletons

items=[]; % J: empty stack
s=1; % J: s := 1;
while s<n % J: while s<n do
    % J: AK = adjacency structure of strong component K with least vertex
    % J: in the subgraph of G induced by {s,s+1,...,n};
    % Use dmperm function to find strong components
    [pp,qq,rr]=dmperm(Ds+eye(n));
    diffr=diff(rr); 
    for aa=1:length(diffr)
        Components(aa,:)=[sort(pp(rr(aa):rr(aa+1)-1)) zeros(1,n-diffr(aa))];
    end
    Components=sortrows(Components,1); % Sort by least vertex
    Components=Components(Components(:,2)>0,:); % Discard singletons
    if isempty(Components)==false % J: if AK not equal to the null set then
        VK=Components(1,:); % Take component K w/ least vertex
        VK=VK(:,VK>0); % Discard zeros
        KeepVK=zeros(n,n);
        KeepVK(VK,VK)=1;
        AK=Ds.*KeepVK;
        s=min(VK); % J: s := least vertex in VK;
        % J: for i in VK do
        blocked(VK)=0; % J: blocked(i) := false;
        B(VK,:)=zeros(length(VK),n); % J: B(i) := null set;        
        dummy=CIRCUIT(s); % J: dummy := CIRCUIT(s);
        Ds(:,1:s)=0; % create subgraph of G induced by {s,s+1,...,n}
        Ds(1:s,:)=0; % create subgraph of G induced by {s,s+1,...,n}
        s=s+1; % J: s := s+1;
    else
        s=n; % J: else s := n;
    end
end
end

function UNBLOCK(u)

global B blocked n

    blocked(u)=0;
    for i=1:n
        B(u,i)=0;
        if blocked(i)==1
            UNBLOCK(i)
        end
    end
end

function f=CIRCUIT(v)

global AK B blocked s items JCycles Ds Fs n

    f=false; % J: f := false;
    items=[items v]; % J: stack v;
    blocked(v)=1; % J: blocked(v) := true;
    
    % *** START MODIFICATION ***
    subcycle=0; 
    if length(items)>2 % Only check if stack could contain a subcycle
        for i=1:(length(items)-2) % For each item
            if Fs(v,items(end-i))==1 % Check if item forms direct subcycle
                subcycle=1;
                break
            end
        end
        if subcycle==0 
            if FloydWarshall(Fs(items,items))==1
                testgroups=nchoosek(items,length(items)-1);
                for i=1:size(testgroups,1)
                    if FloydWarshall(Fs(testgroups(i,:),testgroups(i,:)))==1
                        subcycle=1;
                        break
                    end
                end
            end
        end
    end
     % If a strongly connected component is strictly inside of the stack
    if subcycle==1;
        f=true; 
    else
    % *** END MODIFICATION ***

    AKv=find(AK(v,:)==1);
    for i=1:length(AKv) 
        w=AKv(i); % J: w in AK(v) do
        if w==s  % J: if w=s then
            JCycles{length(JCycles)+1}=[items s]; % J: output circuit
            f=true; % J: f := true;
        elseif blocked(w)==0
            if CIRCUIT(w)==true % J: if CIRCUIT(w)
                f=true; % J: then f := true
            end            
        end
    end
    if f==true
        UNBLOCK(v)
    else
        for i=1:length(AKv)
    	    w=AKv(i); % J: else for w in AK(v) do
            % J: if v not in B(w) then
            B(w,v)=1; % J: put v on B(w)
        end
    end

    end % *** MODIFICATION ***

    items(length(items))=[]; % J: unstack v

end