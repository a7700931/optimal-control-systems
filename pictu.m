function pictu(A,B,R,P,X,TData,G)
u = zeros(size(B,2),length(TData));
Pr = reshape(flipud(P)',size(A,1),size(A,1),[]); % Reshape P into 2 by 2 3D matrix
if nargin <7
    Gr = zeros(length(TData),size(B,1));
else
    Gr = flipud(G);
end
for i=1:length(TData)
    u(:,i) = -(R\B')*(Pr(:,:,i)*X(i,:)'-Gr(i,:)');
end
figure;
plot(TData,u);
title('Optimal Control U(t)'); xlabel('t'); ylabel('U(t)');
legend('U(t)');
end