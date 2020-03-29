# stats 606 project 

# simulation - y (binary time series dataset) x (10 dim time series dataset)

# library package
library(ggplot2)

# scenario1 -----------------------------------------------
# sample size: 1000 with 2 features
# y: binary 
# x1: normal distribution
# x2: random noise 
set.seed(1234)
n = 1000
mu0 = rnorm(1, 0, 1)
mu1 = rnorm(1, 1, 1)
x1 = c(rnorm(n/2, mu0, 1), rnorm(n/2, mu1, 1))
x2 = runif(n)
y = c(rep(0, n/2), rep(1, n/2))
simu_data = data.frame(x1 = x1,  x2 = x2, y = y)

ggplot(simu_data)+
  geom_point(aes(x = x1, y = x2, color = factor(y)))+
  labs(color = 'type')

# bayes error rate p(misclassified) = p(x)p(y neq c(x))dx 
# assume p(y=0) = p(y=1) balanced dataset expand.grid
xnew = t(seq(-3,3,0.01))
x_marginal0 = apply(xnew, 1, function(x) dnorm(x, mu0, 1))
x_marginal1 = apply(xnew, 1, function(x) dnorm(x, mu1, 1))
x_marginal = x_marginal1*0.5 + x_marginal0*0.5
conditional_1 = x_marginal1*0.5/x_marginal
bayes_error = sum(x_marginal*(conditional_1*(conditional_1<0.5) + (1-conditional_1)*(conditional_1>0.5) ))/sum(x_marginal)
bayes_error 


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

dist_scenario1 = distance_cal(simu_data[, 1:2], p = 2, n , t = 0)
dim(dist_scenario1[1,,])
# save and load the dataset
save(simu_data, bayes_error, distance_cal, dist_scenario1, file = 'treeproject.Rdata')

# scenario2 -----------------------------------------------
# sample size: 100 with 10 features
# y: binary 
# x1 - x4: random noise 
# x5 - x10: normal distribution, weights: dirichilet distribution
set.seed(1234)
n = 100
x1 = runif(n)
x2 = runif(n)
x3 = runif(n)
x4 = runif(n)
simu_data = data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
mu0 = rnorm(6, 0, 1)
feature_weight = MCMCpack::rdirichlet(1,rep(1, 6))
for (i in 1:6){
  simu_data[,i+4] = rnorm(n, mu0[i], 1)
}
y_values = as.matrix(simu_data[,5:10], n, 6)%*%t(feature_weight)
threshold = median(y_values)
simu_data$y = (y_values > threshold)*1

ggplot(simu_data)+
  geom_point(aes(x = x1, y = V5, color = factor(y)))+
  labs(color = 'type')

# bayes error rate p(misclassified) = p(x)p(y neq c(x))dx 
# assume p(y=0) = p(y=1) balanced dataset expand.grid
# p(x|y=0) = p(x|w1*x1+...+w10*x10 > 0.5) = p(x1, ..., x10 | z) = p(x1|z)*...p(x10|z)  z also follow normal distribution 
# e.g. p(x = 0 |x+r < 0.5) = p(x=0, r< 0.5)/p(z<0.5)  
# 
xnew = expand.grid(x1=seq(-1,1,0.5), x2=seq(-1,1,0.5), x3=seq(-1,1,0.5), x4=seq(-1,1,0.5), x5=seq(-1,1,0.5), x6=seq(-1,1,0.5))
mean.y = sum(mu0*feature_weight)
var.y = sum(feature_weight^2)
densi0 = function(x){
  x_marginal0 = 1
  for (i in 1:6){
    # same row
    tmp_x_marginal0 = dnorm(x[i], mu0[i], 1)* pnorm(threshold-x[i], mean.y - mu0[i], sqrt(var.y - feature_weight[i]^2))/pnorm(threshold, mean.y, sqrt(var.y))
    x_marginal0 = x_marginal0*tmp_x_marginal0 
  }
  return(x_marginal0)
}
densi1 = function(x){
  x_marginal0 = 1
  for (i in 1:6){
    # same row
    tmp_x_marginal0 = dnorm(x[i], mu0[i], 1)*(1-pnorm(threshold-x[i], mean.y - mu0[i], sqrt(var.y - feature_weight[i]^2)))/(1-pnorm(threshold, mean.y, sqrt(var.y)))
    x_marginal0 = x_marginal0*tmp_x_marginal0 
  }
  return(x_marginal0)
}
x_marginal0 = apply(xnew, 1, densi0)
x_marginal1 = apply(xnew, 1, densi1)
x_marginal = x_marginal1*0.5 + x_marginal0*0.5
conditional_1 = x_marginal1*0.5/x_marginal
bayes_error = sum(x_marginal*(conditional_1*(conditional_1<0.5) + (1-conditional_1)*(conditional_1>0.5) ))/sum(x_marginal)
bayes_error 

dist_scenario = distance_cal(simu_data[, 1:10], p = 10, n , t = 0)
dim(dist_scenario[1,,])
# save and load the dataset
output = list(dist_scenario, simu_data$y)
save(bayes_error, output, file = 'treeproject.Rdata')












# sample size: 1000 with 2 features
# y: binary 
# x1: normal distribution (mu1 = -2, mu2 = 2, sigma = 2)
# x2: random noise 
# gaussian mixture  
# decision boundary 
# bayes optimal error 
set.seed(606)
n = 1000
# mixture of gaussian clusters 
# first we generate 1 means from a bivariate gaussian distribution label it as blue
# then we generate 1 mean from a bivariate gaussian distribution label it as orange 
Sigma = matrix(c(1,0,0,1),2,2)
mu1 = MASS::mvrnorm(n, c(1, 0), Sigma)
mu2 = MASS::mvrnorm(n, c(0, 1), Sigma)
# grid point to decide which is blue or orange 
x1 = c(rnorm(n/2, -1, 1), rnorm(n/2, 1, 1))
x2 = runif(n)
y = c(rep(0, n/2), rep(1, n/2))
simu_data = data.frame(x1 = x1, x2 = x2, y = y)
# bayes error rate 


# scatter plot  x1 vs x2 wrt y 
ggplot(simu_data, aes(x = x1, y= x2, z = y))+
  geom_point(color = simu_data$y)+
  labs(color = 'type')+
  geom_contour(color = 'black')

# tree method illustration
# sample size: 1000 with 2 features
# y: binary 
# x1 & x2: multivariate uniform distribution p =.6
set.seed(606)
n = 1000
#Sigma = matrix(c(10,3,3,2),2,2)
#x = MASS::mvrnorm(n, rep(0, 2), Sigma)
x1 = runif(n)
x2 = runif(n)
# generate y
y = (x1<=0.5)
# overlapping
overlapping = 0.1
x1 = x1 + (x1<=0.5)*overlapping
x1 = x1 - (x1>0.5)*overlapping
simu_data = data.frame(x1 = x1, x2 = x2, y = y)

ggplot(simu_data)+
  geom_point(aes(x = x1, y = x2, color = factor(y)))+
  labs(color = 'type')+
  geom_vline(xintercept = 0.5, color = 'blue')+
  scale_x_continuous(limits = c(0,1))

  

  
# sample size: 1000 with 2 features
# y: binary 
# x1: normal distribution (mu1 = , mu2 = , sigma)
# x2: random noise




# algorithm - find the exemplar, adjust the weight 

# given a time series x matrix (n*p)
# reshape to 3d (ith index, kth dimension, t time series; for each index, t*p dim)
distance_cal = function(x){
  # calculate distance between two exemplar
  # examplar: t*p matrix
  # x: 3d arrary (ith index, kth dimension, t time series; for each index, t*p dim)
  # n: # of time series points we have 
  # w: weight, vector with p*1
  dist = array(rep(0,n*n*p), c(n,n,p))
  
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:p){
        # distance in 3 d: t*p*n  (updated:>>using the huge distance matrix)
        # distance_matrix[i, j, k ] = norm(x[,k,i] - x[,k,j])^2 >> single value to a vector 
        # distance within same dimension (to save time, can only calcualte upper tri-region)
        dist[i,j, k] = norm((x[,k,i] - x[,k,j]), '2')
      }
    }
  }

  # return a vector of index: 1 or 2
  return (dist)
}

classfiy = function(indices_exemplar_1, indices_exemplar_2, w, distance_matrix){
  # calculate distance between two exemplar
  # examplar: t*p matrix
  # x: 3d arrary (ith index, kth dimension, t time series; for each index, t*p dim)
  # n: # of time series points we have 
  # w: weight, vector with p*1
  n = dim(distance_matrix)[1]
  dist = matrix(0L, n, 2)

  for (i in 1:n){
    
    # distance in 3 d: t*p*n  (updated:>>using the huge distance matrix) distance_matrix[i_index,j_index, k]
    # distance_matrix[i, j, k ] = norm(x[,k,i] - x[,k,j])^2 >> single value to a vector 
    dist[i,1] = distance_matrix[indices_exemplar_1, i, ]*w
    dist[i,2] = distance_matrix[indices_exemplar_2, i, ]*w
    
    # minimum distance 
    result[i] = which.min(dist[i,])
  }
  
  # return a vector of index: 1 or 2
  return (result)
}

# return exemaplar and corresponding distance 
# cumulative distance (all other vector to itself with weights) vector, dim: n*1
choose_exemplar = function(region1, distance_matrix, w){
  n_region1 = length(region_1)
  cum_dist = rep(0, n)
  d = length(w)
  
  # cumulative distance (all other vector to itself with weights) vector, dim: n*1
  for (i in 1:n_region1){
    for (j in 1:n_region1){
      for (k in 1:d){
        i_index = region_1[i]
        j_index = region_1[j]
        cum_dist[i] = cum_dist[i] + distance_matrix[i_index,j_index, k]*w[k]
      } 
    }
  }
  
  # return exemplar index (original order) and cum_dist
  index_exemplar = which.min(cum_dist)
  index_exemplar = region_1[index_exemplar]
  #distance_stored = cum_dist[index_exemplar]
  
  #return(c(index_exemplar, distance_stored))
  return(index_exemplar)
}

# find the best weight 
# return a vector 
find_weight = function(indices_exemplar_1, misclassified_1, distance_matrix, w){
  p = length(w)
  mis_dis = rep(0, p)
  i_index = indices_exemplar_1
  n_region1_mis = length(misclassified_1) 
  
  for (k in 1:p){
    for (j in 1:n_region1_mis){
      j_index = misclassified_1[j]
      mis_dist[k] = mis_dist[k] + distance_matrix[i_index,j_index, k]
    } 
    # update - comparing with 0
    mis_dist[k] = max(mis_dist[k], 0)
  }
  # normizaltion - return a vector 
  w =  mis_dist/sum(mis_dist)
  
  return(w)
}

exemplar_weight = function(indices_exemplar_1, indices_exemplar_2, w, distance_matrix){
  
  # classify the rest points into two areas - based on distance (lie within certain distance)
  # go through everypoint into the area to classify points into 2 area 
  index_classification = classfiy(indices_exemplar_1, indices_exemplar_2, w, distance_matrix)
  
  # find the exemplar - same response with minimize distance
  # two region: 
  # in region 1, index who were classfied correctly, return a vector of index 
  region_1 = which(index_classification == 1 & y == 0)
  region_2 = which(index_classification == 2 & y == 1)
  
  # calculate distance array (same type) between x values  3d: i:row index, j: column column, k: dimension
  # we have distance matrix: stored all pointwise distance between 2 points from same group
  # distance_matrix_1 = array(distance_data, c(n, n, d))  x[,,i] t*p dimension for index i 
  # distance_matrix[i, j, k ] = norm(x[,k,i] - x[,k,j])^2 >> single value to a vector 

  # update exemplar
  indices_exemplar_1 = choose_exemplar(region1, distance_matrix, w)
  indices_exemplar_2 = choose_exemplar(region2, distance_matrix, w)
  exemplar_1 = x[,,indices_exemplar_1]
  exemplar_2 = x[,,indices_exemplar_2]
  
  # update the weight for region 1 
  # classfication with new exemplar
  index_classification = classfiy(indices_exemplar_1, indices_exemplar_2, w, distance_matrix)
  
  # given exemplar and distance matrix, maxmize distance to whom are misclassfied 
  misclassified_1 = which(index_classification == 1 & y == 1)
  misclassified_2 = which(index_classification == 2 & y == 0)
  
  # algorthms to find the best weight for each k 
  w_1 =  find_weight(indices_exemplar_1, misclassified_1, distance_matrix, w)
  w_2 =  find_weight(indices_exemplar_2, misclassified_2, distance_matrix, w)
  
  # return index for exemplar and respective weight vector
  return(c(indices_exemplar_1, indices_exemplar_2,  w_1, w_2))
}



# initization 
# x 3d array
# x is 3d arrary [,,i] ith index -- t*p time series 
t = dim(x)[1]
p = dim(x)[2]
n = dim(x)[3]

# given a distance_matrix
distance_matrix

# weight is vector dim: 1*k
w = rep(1/d, p)
  
# two groups in y
y_group0 = which(y==0)
y_group1 = which(y==1)

# random two exemplar, classify datapoints into two areas (maybe chosen exemplar with same type) >> chosen from two different index set
indices_exemplar_1 = sample(y_group0, 2)
indices_exemplar_2 = sample(y_group1, 2)

exemplar_1 = x[,,indices_exemplar[1]]
exemplar_2 = x[,,indices_exemplar[2]]

# for-loop to find the best exemplar and weight 
#?? tree structure, each time return a different vector of weight

results = exemplar_weight(exemplar_1, exemplar_2, w, distance_matrix)
indices_exemplar_1 = results[1]
indices_exemplar_2  = results[2]
w_1 = results[3]
w_2 = results[4]
# stopping criterion on weight, see how the difference in weight change sum(w1-w0)/sum(w1)
# weight can be negative or positive 








