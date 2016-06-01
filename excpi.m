N = 50

%pi(x) = boxcar function
x = -pi:pi
p = zeros(size(x))

region1 = x(x < -pi/2 | x > pi /2)
region2 = x(x >= -pi/2 & x <= pi/2)



for i = 0:1:size(x)
    for n = -N:N
        %compute partial sum
        p[i] = p[i] + 
   
    end
end