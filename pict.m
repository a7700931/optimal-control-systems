function pict(ICsize,TData,Data,z)
num = 1:numel(ICsize);
index = reshape(num,size(ICsize))';
Data(:,tril(index,-1)~=0)=[];
[row,col] = find(triu(index)==index);
figure,plot(TData,Data);
xlabel('t'); ylabel(sprintf('%s(t)',inputname(3)));

if 'P'== inputname(3)
    title('Riccati coefficients P(t)');
    legend(compose("P%d%d(t)",row,col));
elseif nargin == 3 && 'X' == inputname(3)
    title('Optimal States X(t)');
    legend(compose("X%d(t)",col));
elseif nargin == 4 && 'X' == inputname(3)
    title('State X(t)');
    hold on 
    plot(TData,z(1,:),'--k',TData,z(2,:),':k','LineWidth',0.5)
    lgd = legend([compose("X%d(t)",col),'Z1(t)','Z2(t)']);
    lgd.NumColumns = 2;
    hold off
else
    title('Vector G(t)');
    legend(compose("G%d(t)",col));    
end
end