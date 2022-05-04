dist01 = randn(100,1)*0.5 + 4;
dist02 = randn(100,1)*0.5 + 4.8;

cat01 = zeros( length(dist01), 1 );
cat02 = zeros( length(dist02), 1 ) + 1;

pred = [dist01; dist02 ]; % concatenate the two groups
resp = [cat01; cat02 ]; % create the dichotomous variable for the two groups (1 and 0)

mdl = fitglm( pred, resp, 'Distribution', 'binomial', 'Link', 'logit' );% fit
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(resp,scores,1); %roc


figure(1)
subplot(1,3,1); hist(dist01);
subplot(1,3,2); hist(dist02);
subplot(1,3,3); 
AUC
plot(X,Y)
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')
