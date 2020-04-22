function err=error_vectors(y,yhat)
ReMSE = mean(abs(y - yhat))/max(y);  % Root Mean Squared Error

RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

err=[ReMSE, RMSE];
end