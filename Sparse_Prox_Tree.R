#######################

## Sparse Proximity Tree

#######################

## Overview of the algorithm
## Given (X_i, Y_i), i = 1,2,...,n as the data points, where X_i could be high-dimensional vectors/time-series
## Y_i as a binary label that can only take value of 0 or 1, our algorithm is aiming at classifying the data points.
## Our algorithm will choose an exemplar from each of the two classes. All samples will be assigned to one of the 
## exemplars depending on which exemplar is closer to the sample. We will update the choice of exemplars to have
## the best classfication result given a distance metric. Then we fix the exemplar, and update the distance metric.
## The final output of our algorithm are two exemplars and a weight vector for all features at a decision node.

library(data.tree)
load("treeproject.Rdata")

# Function: create_dobj
#' @details create a data object for sparse proximity tree
#' @param data: a list containing the pairwise distance matrix and label
#' @param iter: number of iterations for (clustering)-(soft-thresholding) cycle
#' @param s: tuning parameter for weight vector
#' @param index: indices of sampels that are waiting to be classified at the node

create_dobj = function(data,iter,s,index,layer,converge = 0){
  dist_data = data$dist_scenario # distance matrix
  y = data$y # cluster label
  p = dim(dist_data)[1] # number of features
  N = dim(dist_data)[2] # number of samples
  w = rep(1/sqrt(p),p) # initialize the weight parameter
  #index = 1:N # indices for all samples waiting to be classfied at the node
  
  spar_prox = list(p = p, N = N, iter = iter, s = s, index = index, w = w, layer = layer, converge = converge) # data object for sparse proximity tree
  return(spar_prox)
}

# Function: dist_cal
#' @details: this function create a weighted distance matrix
#' @param: weight: feature weight
#' @param: dist_data: feature-wise distance information

dist_cal = function(weight,dist_data){
  p = length(weight)
  N = dim(dist_data)[2]
  dist_mat = matrix(0,nrow=N,ncol=N)
  
  for (i in 1:p){
    dist_mat = dist_mat + weight[i]*matrix(dist_data[i,,],nrow=N,ncol=N)
  }
  return(dist_mat)
}

# Function: soft_thres
#' @details: this function does soft-thresholding for vector vec given l1-nrom constraint s and l2-norm constraint 1
#' @param: vec: vector which will be soft-thresholded
#' @param: s: tuning parameter for penalty l1-norm of feature weight
#' @param: delta: constant for adjusting weights for feasibility
#' @param: lr: learning rate of delta

soft_thres = function(vec,s,delta=0,lr=0.01){
  vec = as.matrix(vec,ncol=1)
  a_plus = as.matrix(apply(vec,1,function(x) max(x,0)),ncol=1)
  S = sign(a_plus)*as.matrix(apply(a_plus,1,function(x) max(abs(x)-delta,0)),ncol=1)
  w = S/norm(S,"2")
  
  
  if (norm(w,"1")>s){
    while(FALSE){
      delta = delta + lr
      S = sign(a_plus)*as.matrix(apply(a_plus,1,function(x) max(abs(x)-delta,0)),ncol=1)
      if (norm(S,"1")<=s*norm(S,"2")){
        w = S/norm(S,"2")
        break
      }
    }
  }
 
  return(w) 
}

# Function: node_split
#' @details Given data, iteration steps and tuning parameter, select two exemplars from the data and attain weights
#' @param data: a list containing the pairwise distance matrix and label
#' @param spar_prxo: node specific information

node_split = function(data,spar_prox,lr=0.01,epsilon=1e-4){
  ind = spar_prox$index # find samples that await classification at the node
  #dist_data = data$dist_scenario[,ind,ind] # retrieve the distance information of the samples relevant
  
  y0 = intersect(which(data$y==0),ind) # find the indices for all samples of class 0 that are alive
  y1 = intersect(which(data$y==1),ind) # find the indices for all samples of class 1 that are alive
  w = spar_prox$w # feature weight parameter
  s = spar_prox$s # tuning parameter
  p = spar_prox$p # number of features
  
  iter = spar_prox$iter # number of iterations
  iter_counter = 1
  
  converge = spar_prox$converge # the convergence status from the previous node
  exitflag = F
  
  while(iter_counter<=iter & exitflag == F){
    # initialize the medoid for class 0 and class 1
    dist_mat = dist_cal(weight = w, dist_data = data$dist_scenario) # distance matrix of all samples
    
    
    if (length(y0)==1){
      medoid_0 = y0
    }else{
      dist0 = dist_mat[y0,y0]
      r0 = rowSums(dist0)
      v0 = colSums(apply(dist0, 2, function(x) x/r0))
      medoid_0 = min(y0[which(v0==min(v0))]) # initial medoid for cluster 0
    }
    
    if (length(y1)==1){
      medoid_1 = y1
    }else{
      dist1 = dist_mat[y1,y1]
      r1 = rowSums(dist1)
      v1 = colSums(apply(dist1, 2, function(x) x/r1))
      medoid_1 = min(y1[which(v1==min(v1))]) # initial medoid for cluster 1
    }
    
    # now we are going to start the internal loop of finding the best medoid under w
    for (i in 1:iter){
      # assign points to two medoids
      cluster0 = intersect(which(dist_mat[,medoid_0]<dist_mat[,medoid_1]),ind)
      cluster1 = intersect(which(dist_mat[,medoid_0]>=dist_mat[,medoid_1]),ind)
      
      # update the medoid
      clust_y0 = intersect(cluster0,y0)
      if (length(clust_y0)==1){
        new_medoid_0 = medoid_0
      }else{
        dist_sum0 = rowSums(dist_mat[clust_y0,clust_y0])
        new_medoid_0 = min(clust_y0[which(dist_sum0==min(dist_sum0))])
      }
      
      
      clust_y1 = intersect(cluster1,y1)
      if (length(clust_y1)==1){
        new_medoid_1 = medoid_1
      }else{
        dist_sum1 = rowSums(dist_mat[clust_y1,clust_y1])
        new_medoid_1 = min(clust_y1[which(dist_sum1==min(dist_sum1))])
      }
      
      
      if (medoid_0==new_medoid_0 & medoid_1==new_medoid_1){
        break
      }else{
        medoid_0 = new_medoid_0
        medoid_1 = new_medoid_1
      }
    }
    
    # now we update the weights using soft-thresholding
    # the first thing is to calculate a for each feature
    a = NULL
    cluster0 = intersect(which(dist_mat[,medoid_0]<dist_mat[,medoid_1]),ind)
    cluster1 = intersect(which(dist_mat[,medoid_0]>=dist_mat[,medoid_1]),ind)
    cluster_correct = c(intersect(cluster0,y0),intersect(cluster1,y1))
    
    for (j in 1:p){
      delta_d = data$dist_scenario[j,,medoid_0] - data$dist_scenario[j,,medoid_1]
      delta_y = ifelse(data$y==1,1,-1)
      delta_y[!cluster_correct] = 0
      a = c(a,crossprod(delta_d[ind],delta_y[ind]))
    }
    
    # check if a has negative sign for all elements, we break and do not split anymore
    if (max(a)<=0){
      converge = 2
      break
    }
    
    # update weight
    w = soft_thres(vec = a,s = s)
    
    # now we check if convergence criteria are fulfilled:
    # 1. medoids are no longer changing
    # 2. weights are changing very little
    
    if (iter_counter==1){
      old_medoid_0 = medoid_0
      old_medoid_1 = medoid_1
      old_w = w
    }else{
      if (old_medoid_0 == medoid_0 & old_medoid_1 == medoid_1 & norm(old_w-w,"2")<epsilon){
        exitflag = T
      }else{
        old_medoid_0 = medoid_0
        old_medoid_1 = medoid_1
        old_w = w
      }
    }
    iter_counter = iter_counter + 1
  }
  
  # now we shall create two new spar_prox objects and also return the medoid and weight information
  dist_mat = dist_cal(weight = w, dist_data = data$dist_scenario) # distance matrix of all samples
  cluster0 = intersect(which(dist_mat[,medoid_0]<dist_mat[,medoid_1]),ind)
  cluster1 = intersect(which(dist_mat[,medoid_0]>=dist_mat[,medoid_1]),ind)
  
  new_converge = ifelse(converge!=2,0,2)
  current_layer = spar_prox$layer
  node0 = create_dobj(data = data,iter = iter,s = s,index = cluster0, layer = current_layer+1,converge = new_converge)
  node1 = create_dobj(data = data,iter = iter,s = s,index = cluster1, layer = current_layer+1,converge = new_converge)
  node_info = list(medoid_0 = medoid_0, medoid_1 = medoid_1, w = w, iter_step = iter_counter)
  
  return(list(node0 = node0, node1 = node1, node_info = node_info))
}

# Function: node_create
#' @details This is a functional wrapper of the node_split that creates a root/branch/leaf node of the sparse proximity tree
#' @param parent: the parent node
#' @param other parameters follow the definitions in node_split

node_create = function(parent, data, spar_prox, lr=0.01,epsilon=1e-4){
  # check if the spar_prox object is pure
  y = data$y
  ind = spar_prox$index
  
  if (length(unique(y[ind])) == 1){
    child = parent$AddChild(paste("class",unique(y[ind]))) # create a child node if the node is pure
    
    # if pure, write the node information into the leaf node features
    child$class0_count = length(which(y[ind] == 0)) 
    child$class1_count = length(which(y[ind] == 1))
    child$obs_count = child$class0_count + child$class1_count
    child$Gini_index = 1 - ((child$class0_count)^2 + (child$class1_count)^2)/(child$obs_count)^2 
    parent$layer = spar_prox$layer
    child$layer = spar_prox$layer+1
    child$converge = 0
    parent$converge = spar_prox$converge
    parent$majority = ifelse(length(which(y[ind] == 0))>length(which(y[ind] == 1)),0,1)
    
  }else{
    # if not pure, continue on splitting
    result = node_split(data = data, spar_prox = spar_prox, lr = lr, epsilon = epsilon)
    
    # write node0 and node 1 into child node
    child0 = parent$AddChild(result$node_info$medoid_0)
    child1 = parent$AddChild(result$node_info$medoid_1)
    parent$weight = result$node_info$w
    parent$majority = ifelse(length(which(y[ind] == 0))>length(which(y[ind] == 1)),0,1)
    parent$layer = spar_prox$layer
    parent$converge = result$node0$converge
      
    node_create(parent = child0, data = data, spar_prox = result$node0, lr = lr, epsilon =  epsilon)
    node_create(parent = child1, data = data, spar_prox = result$node1, lr = lr, epsilon =  epsilon)
  }
}


# Function: build_tree
#' @details this function returns a sparse proximity tree classifier
#' @param data: a list containing the pairwise distance matrix and label
#' @param lr: step length of the adjustment of delta in soft-thresholding
#' @param epsilon: convergence criterion for weight estimator
#' @param depth: maximum depth of the sparse proximity tree

build_tree = function(data,iter = 100, lr=0.01,epsilon=1e-4, depth=3, train_ind = NULL){
  N = dim(data$dist_scenario)[2]
  p = dim(data$dist_scenario)[1]
  if (is.null(train_ind)){
    spar_prox = create_dobj(data = data, iter = iter, s = sqrt(p), index = 1:N, layer = 0, converge = 0)
  }else{
    spar_prox = create_dobj(data = data, iter = iter, s = sqrt(p), index = train_ind, layer = 0, converge = 0)
  }
  
  tree = Node$new("Sparse Proximity Tree")
  node_create(parent = tree, data = data, spar_prox = spar_prox, lr = lr, epsilon = epsilon)
  treeClone = Clone(tree, pruneFun = function(x) (x$layer<=depth))

  
  treeClone$Do(function(x) if(!isNotLeaf(x)) {x$AddChild(paste("Class", x$majority))})
  treeClone$Do(function(x) if(!isNotLeaf(x)) {x$converge=0})
  return(treeClone)
}


# Function: tree_test
#' @details this function gives prediction on a test sample with a fitted tree
#' @param data: training data information (not the distance matrix, but the true features)
#' @param tree: fitted sparse proximity tree object
#' @param feature: test sample input feature

tree_test = function(data,tree,feature,dist_func = NULL){
  # if it reaches a leaf node, stop and output the majority class at the leaf node
  if (tree$leafCount==1){return(tree$majority)} 
  
  # otherwise we have to decide where to go
  w = tree$weight
  medoid_0 = as.integer(tree$children[[1]]$name)
  medoid_1 = as.integer(tree$children[[2]]$name)
  p = ncol(data) - 1
  x0 = data[medoid_0,1:p]
  x1 = data[medoid_1,1:p]
  feature = matrix(feature,ncol = 1)
  
  if (is.null(dist_func)){
    # use Euclidean distance
    d0 = norm(x0-feature,"2")
    d1 = norm(x1-feature,"2")
  }else{
    d0 = dist_func(x0,feature)
    d1 = dist_func(x1,feature)
  }
  
  if (d0<d1){
    return(tree_test(data = data, tree = tree$children[[1]], feature = feature, dist_func = dist_func))
  }else{
    return(tree_test(data = data, tree = tree$children[[2]], feature = feature, dist_func = dist_func))
  }
}

# build tree for a simulation scenario
# train test split
input_data = output1
data = list(dist_scenario = input_data[[1]], y = input_data[[2]])
tree = build_tree(data = data,depth = 3,train_ind = train_ind,iter = 99)
plot(tree)

# Function: tree_eval
#' @details this function does train test split on the input data and provide prediction accuracy on the test set using Sparse Proximity Tree
#' @param input_data: a list of three components: distance matrix (p by N by N), binary label, feature-label dataset
#' @param train_size: proportion of samples used for training the sparse proximity tree

tree_eval = function(input_data,train_size = 0.7){
  input_data = output1
  N = dim(input_data[[1]])[2]
  p = dim(input_data[[1]])[1]
  train_ind = sort(sample(1:N,size = floor(train_size*N),replace = F))
  test_ind = sort(setdiff(1:N,train_ind))
  
  data = list(dist_scenario = input_data[[1]], y = input_data[[2]])
  tree = build_tree(data = data,depth = 3,train_ind = train_ind)
  
  # test it!
  D = input_data[[3]]
  test_pred = apply(D[test_ind,1:p],MARGIN = 1,function(x) tree_test(data = D, tree = tree, feature = x))
  acc = sum(test_pred==data$y[test_ind])/length(test_ind)
  
  return(acc)
}


alldata = list(output1,output2,output3,output4,output4_sd)
acc = lapply(alldata,tree_eval)


