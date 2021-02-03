function rmse=calc_rmse(a,b)

rmse=sqrt(mean((a(:)-b(:)).^2))*100;
