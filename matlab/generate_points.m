%

clear;





n_points = [10 10 10 0];
dx       = 1./n_points; 

num = 1;
for i=0:n_points(1)+1
    for j=0:n_points(2)-i
        for k=0:n_points(3)-j-i
            l = 10-i-j-k;
            
            prop(num,:) = [dx(1)*i dx(2)*j dx(3)*k dx(4)*l];
            
            num=num+1;
        end
    end
end
            