x1 = (0:0.05:1)';
x2 = (0:0.05:1)';

nx1 = length(x1);
nx2 = length(x1);

knot = [0.25 0.25; 0.25 0.5;0.75 0.25;0.75 0.75];
thinPlate = zeros(nx1,nx2,size(knot,1));
hardy = zeros(nx1,nx2,size(knot,1));

for k = 1:size(knot,1)
    for i = 1:nx1
        for j = 1:nx2
            x = [x1(i) x2(j)]';
            u = norm(x-knot(k,:)');
            if u>0
                thinPlate(i,j,k) = log(u)*(u^2);
            end
            hardy(i,j,k) = sqrt(0.1^2 + u^2);
        end
    end
end

[X Y] = meshgrid(x1,x2);


figure
hold on
%subplot(1,2,1)
for k = 1:size(knot,1)
    %subplot(2,2,k)
    surf(X,Y,sum(thinPlate(:,:,:),3))
    %title(['knot = [',num2str(knot(k,1),2),' ',num2str(knot(k,2),2),']'])
end
%subplot(1,2,2)
%surf(X,Y,hardy)
