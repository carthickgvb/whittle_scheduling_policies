function [out out1] = findMIS_new(G)

n = numnodes(G);
temp = de2bi(2^n - 1);
A = adjacency(G);
out1=zeros(2^n - 1,n);%%%%%%
out = 0;
for i=1:2^n -1
    idx = logical(de2bi(i,size(temp,2)));
    a = A(idx,idx);
    if(a == zeros(size(a)) & size(a,1) >= out)
        out = size(a,1);
        out1(i,:) =logical(idx);%%%%%%%
    end
end
out1( all(~out1,2), : ) = [];
out1=flipud(out1);
rw=size(out1,1);
for k=1:rw-1
  for h=k+1:rw
     if(~any(out1(k,:)-out1(h,:)<0))
      out1(h,:)=0; 
     end
  end
end
out1( all(~out1,2), : ) = [];
end