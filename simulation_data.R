# stats 606 project 

# simulation - y (binary time series dataset) x (10 dim time series dataset)

# library package
library(ggplot2)

# scenario1 -----------------------------------------------
# sample size: 1000 with 2 features
# y: binary 
# x1: normal distribution
# x2: random noise 
set.seed(123)
n = 200
mu0 = rnorm(1, 0, 1)
mu1 = rnorm(1, 1, 1)
x1 = c(rnorm(n/2, mu0, 1), rnorm(n/2, mu1, 1))
x2 = runif(n)
y = c(rep(0, n/2), rep(1, n/2))
simu_data = data.frame(x1 = x1,  x2 = x2, y = y)

ggplot(simu_data)+
  geom_point(aes(x = x1, y = x2, color = factor(y)))+
  labs(color = 'class')

# bayes error rate p(misclassified) = p(x)p(y neq c(x))dx 
# assume p(y=0) = p(y=1) balanced dataset expand.grid
xnew = t(seq(-3,3,0.01))
x_marginal0 = apply(xnew, 1, function(x) dnorm(x, mu0, 1))
x_marginal1 = apply(xnew, 1, function(x) dnorm(x, mu1, 1))
x_marginal = x_marginal1*0.5 + x_marginal0*0.5
conditional_1 = x_marginal1*0.5/x_marginal
bayes_error1 = sum(x_marginal*(conditional_1*(conditional_1<0.5) + (1-conditional_1)*(conditional_1>0.5) ))/sum(x_marginal)
bayes_error1 


# distance matrix p*n*n
distance_cal = function(xdata, p, n, t){
  # x: 3d arrary (ith index, kth dimension, t time series; for each index, t*p dim)
  # n: # of time series points we have 
  m = n-t 
  dist = array(rep(0,p*n*n), c(p,n,n))
  x = array(rep(0,p*t*n), c(p,t+1,n))
  for (i in 1:n){
    for (k in 1:p){
      # time series dataset x[k, ,i]: vector 1:t 
      x[k, ,i] = xdata[i:i+t,k]
    }
  }
  
  for (i in 1:m){
    for (j in 1:m){
      for (k in 1:p){
        # x in 3 d: t*p*n  (updated:>>using the huge distance matrix)
        # distance_matrix[i, j, k ] = norm(x[,k,i] - x[,k,j])^2 >> single value to a vector 
        # distance within same dimension (to save time, can only calcualte upper tri-region)
        dist[k,i,j] = norm((x[k, ,i] - x[k, ,j]), '2')
      }
    }
  }
  return (dist)
}

dist_scenario = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
# save and load the dataset
output1 = list(dist_scenario, simu_data$y, simu_data)
#save(simu_data, bayes_error, distance_cal, dist_scenario1, file = 'treeproject.Rdata')


# scenario3 (consider heteroskedsity inside same group) -----------------------------------------------
# sample size: 1000 with 2 features
# y: binary 
# x1: normal distribution
# x2: random noise 
n = 200
x1 = runif(n)

x2 = c(lapply(x1[1:(n/4)], function(x) runif(1, -0.15, x-0.15) )%>%unlist(),
       lapply(x1[(n/4+1):(n/2)], function(x) runif(1, x+0.15, 1.15))%>%unlist(),
       lapply(x1[(n/2+1):n], function(x) runif(1, x-0.15, x+0.15))%>%unlist())
              
y = c(rep(0, n/2), rep(1, n/2))
simu_data = data.frame(x1 = x1,  x2 = x2, y = y)

ggplot(simu_data)+
  geom_point(aes(x = x1, y = x2, color = factor(y)))+
  labs(color = 'class')

# bayes error rate p(misclassified) = p(x)p(y neq c(x))dx 
# assume p(y=0) = p(y=1) balanced dataset expand.grid 
# xnew = t(seq(-3,3,0.01))
# xnew = t(simu_data[,1])
# x_marginal0 = apply(xnew, 1, function(x) dnorm(x, mu0, 1)*(dnorm(x, mu0, 1)>dnorm(x, mu2,1)) + dnorm(x, mu2,1)*(dnorm(x, mu0, 1)<dnorm(x, mu2,1)))
# x_marginal1 = apply(xnew, 1, function(x) dnorm(x, mu1, 1))
# x_marginal = x_marginal1*0.5 + x_marginal0*0.5
# conditional_1 = x_marginal1*0.5/x_marginal
# bayes_error2 = sum(x_marginal*(conditional_1*(conditional_1<0.5) + (1-conditional_1)*(conditional_1>0.5) ))/sum(x_marginal)
# bayes_error2 

dist_scenario = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
# save and load the dataset
output3 = list(dist_scenario, simu_data$y, simu_data)


# scenario2 -----------------------------------------------
# sample size: 100 with 10 features
# y: binary 
# x1 - x4: random noise 
# x5 - x10: normal distribution, weights: dirichilet distribution
set.seed(123)
n = 200
mu0 = rnorm(1, 0, 1)
mu1 = rnorm(1, 1, 1)
x1 = c(rnorm(n/2, mu0, 1), rnorm(n/2, mu1, 1))
x2 = runif(n)
y = c(rep(0, n/2), rep(1, n/2))
simu_data = data.frame(x1 = x2,  x2 = x1, y = y)

# rotation matrix 
theta1 = -pi/8
theta2 = pi/4
r = function(theta ) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)

for(i in 1:n){
  if(i < (n/2)){
    simu_data[i,1:2] = as.matrix(simu_data[i,1:2], 1, 2)%*%r(theta1)
  }
  else{
    simu_data[i,1:2] = as.matrix(simu_data[i,1:2], 1, 2)%*%r(theta2)
  }
}

ggplot(simu_data)+
  geom_point(aes(x = x1, y = V5, color = factor(y)))+
  labs(color = 'class')

# bayes error rate p(misclassified) = p(x)p(y neq c(x))dx 
# assume p(y=0) = p(y=1) balanced dataset expand.grid
# p(x|y=0) = p(x|w1*x1+...+w10*x10 > 0.5) = p(x1, ..., x10 | z) = p(x1|z)*...p(x10|z)  z also follow normal distribution 
# e.g. p(x = 0 |x+r < 0.5) = p(x=0, r< 0.5)/p(z<0.5)  
# 
dist_scenario = distance_cal(simu_data[, 1:2], p = 10, n , t = 0)
# save and load the dataset
output2 = list(dist_scenario, simu_data$y, simu_data)




# scenario4 -----------------------------------------------
# sample size: 100 with 10 features
# y: binary 
# x1 - x2: random noise (uniform)
# x3: random noise (log-normal)
# x4: heavy tail distribution: log-normal
# x5 - x10: normal distribution, 
# weights: dirichilet distribution
set.seed(123)
n = 200
x1 = runif(n)
x2 = runif(n)
simu_data = data.frame(x1 = x1, x2 = x2)
mu0 = rnorm(8, 0, 1)
feature_weight = MCMCpack::rdirichlet(1,rep(1, 7))
for (i in 1:8){
  if (i<3){
    simu_data[,i+2] = rlnorm(n, mu0[i], 1)
  }
  else{
    simu_data[,i+2] = rnorm(n, mu0[i], 1) 
  }
}
y_values = as.matrix(simu_data[,4:10], n, 7)%*%t(feature_weight)
threshold = median(y_values)
simu_data$y = (y_values > threshold)*1

ggplot(simu_data)+
  geom_point(aes(x = x1, y = V4, color = factor(y)))+
  labs(color = 'class')



# without standizarion
dist_scenario = distance_cal(simu_data[, 1:10], p = 10, n , t = 0)
dim(dist_scenario[1,,])
# save and load the dataset
output4 = list(dist_scenario, simu_data$y, simu_data)

# with standizarion (train , test split!!!)
index = sample(1:200, 30)
train = simu_data[-index, ]
test = simu_data[index, ]
mean_train = colMeans(train[,-11])
std_train = sqrt(diag(var(train[,-11])))
train_x = scale(train[,-11], center = mean_train, std_train)
all_x = scale(simu_data[,-11], center = mean_train, std_train)
dist_scenario_train = distance_cal(train_x, p = 10, n-30 , t = 0)
dist_scenario_all = distance_cal(all_x, p = 10, n , t = 0)
# save and load the dataset
output4_sd = list(dist_scenario_train, train$y, dist_scenario_all, train, test)
save(bayes_error1, output1, bayes_error2, output2, bayes_error3, output3, 
     bayes_error4, output4, output4_sd, file = 'treeproject.Rdata')



# 5*2 cv test
# corrected repeated 5-fold cross validation test 

scenario1 = function(){
  #set.seed(123)
  n = 200
  mu0 = rnorm(1, 0, 1)
  mu1 = rnorm(1, 1, 1)
  x1 = c(rnorm(n/2, mu0, 1), rnorm(n/2, mu1, 1))
  x2 = runif(n)
  y = c(rep(0, n/2), rep(1, n/2))
  simu_data = data.frame(x1 = x1,  x2 = x2, y = y)
  # random shuffle the points 
  index = sample(1:n, n)
  simu_data = simu_data[index,]
  dist_scenario = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
  output1 = list(dist_scenario, simu_data$y, simu_data)
  return(output1)
}  


scenario2 = function(){
  #set.seed(123)
  n = 200
  mu0 = rnorm(1, 0, 1)
  mu1 = rnorm(1, 1, 1)
  x1 = c(rnorm(n/2, mu0, 1), rnorm(n/2, mu1, 1))
  x2 = runif(n)
  y = c(rep(0, n/2), rep(1, n/2))
  simu_data = data.frame(x1 = x1,  x2 = x2, y = y)
  
  # rotation matrix 
  theta1 = -pi/8
  theta2 = pi/4
  r = function(theta ) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)
  
  for(i in 1:n){
    if(i < (n/2)){
      simu_data[i,1:2] = as.matrix(simu_data[i,1:2], 1, 2)%*%r(theta1)
    }
    else{
      simu_data[i,1:2] = as.matrix(simu_data[i,1:2], 1, 2)%*%r(theta2)
    }
  }
  # random shuffle the points 
  index = sample(1:n, n)
  simu_data = simu_data[index,]
  dist_scenario = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
  output2 = list(dist_scenario, simu_data$y, simu_data)
  
  return(output2)
}  

scenario3 = function(){
  n = 200
  x1 = runif(n)
  x2 = c(lapply(x1[1:(n/4)], function(x) runif(1, -0.15, x-0.15) )%>%unlist(),
         lapply(x1[(n/4+1):(n/2)], function(x) runif(1, x+0.15, 1.15))%>%unlist(),
         lapply(x1[(n/2+1):n], function(x) runif(1, x-0.15, x+0.15))%>%unlist())
  
  y = c(rep(0, n/2), rep(1, n/2))
  simu_data = data.frame(x1 = x1,  x2 = x2, y = y)
  # random shuffle the points 
  index = sample(1:n, n)
  simu_data = simu_data[index,]
  dist_scenario = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
  output3 = list(dist_scenario, simu_data$y, simu_data)
  
  return(output3)
}  
                
                
result_eval = function(input_data, K = 5){
    N = dim(input_data[[1]])[2]
    p = dim(input_data[[1]])[1]
    
    # 5-fold cross validation
    fold_size = floor(N/K)
    acc = rep(0, K)
    acc2 = rep(0, K)
    for( i in 1:K){
      if (i!= K){
        cv_ind = ((i-1)*fold_size+1) : (i*fold_size)
      }
      else{
        cv_ind = ((i-1)*fold_size+1):N
      }
      
      train_ind = sort(setdiff(1:N,cv_ind))
    
      data = list(dist_scenario = input_data[[1]], y = input_data[[2]])
      tree = build_tree(data = data,depth = 4,train_ind = train_ind)
      # cart tree
      cartt.model = rpart(y~., data = input_data[[3]][train_ind, ], method = 'class', maxdepth = 3)
      
      # test it!
      D = input_data[[3]]
      test_pred = apply(D[cv_ind,1:p],MARGIN = 1,function(x) tree_test(data = D, tree = tree, feature = x))
      test_pred2 = predict(cartt.model, input_data[[3]][cv_ind, ], type = 'class')
      acc[i] = sum(test_pred==data$y[cv_ind])/length(cv_ind)
      acc2[i]= sum(test_pred2==data$y[cv_ind])/length(cv_ind)
    }

    return(list(acc, acc2, acc-acc2))
}

                        
k = 5
r = 100
n = 200
                        

results = data.frame(matrix(rep(0, k*r*6), k*r, 6))
names(results) = c('x1', 'x2', 'y1', 'y2', 'z1', 'z2')
acc.dif = data.frame(matrix(rep(0, k*r*3), k*r, 3))

set.seed(123)
for (i in 1:r){
  # average accuracy rate 
  alldata = list(scenario1(),scenario2(),scenario3())
  acc = lapply(alldata,result_eval)
  # accuracy difference
  for (j in 1:3){
    results[(5*i-4):(5*i),2*j-1] = acc[[j]][[1]]     
    results[(5*i-4):(5*i),2*j] = acc[[j]][[2]]
    acc.dif[(5*i-4):(5*i),j] = acc[[j]][[3]]
  }
}

m = colMeans(acc.dif)
sigma2 = apply(acc.dif, 2, var)
t = m/sqrt((1/k/r+0.2/0.8)*sigma2)
t   #1.96                        
