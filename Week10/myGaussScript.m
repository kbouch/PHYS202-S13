x = -10:0.1:5;
sigma = 2;
y = exp(-x.^2/sigma);
plot(x,y)
xlabel('x-axis (m)')
ylabel('Gaussian')
title('myGaussian')