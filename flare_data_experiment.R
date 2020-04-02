# apply proximity tree to Solar Flare data
library(reticulate)
np = import("numpy")

D = np$load("Spar_dist.npy")
y_train = np$load("train_label.npy")
y_test = np$load("test_label.npy")
y = c(y_train,y_test)
y[which(y==2)] = 1

feature_name = c('TOTUSJH', 'log(TOTPOT)', 'TOTUSJZ', 'ABSNJZH', 'log(SAVNCPP)', 'log(USFLUX)',
                 'log(AREA_ACR)', 'MEANPOT', 'R_VALUE', 'SHRGT45', 'MEANSHR', 'MEANGAM', 'MEANGBT', 'MEANGBH',
                 'MEANGBZ', 'MEANJZH', 'MEANJZD', 'MEANALP')



train_samples = 1:length(y_train)
flare_data = list(dist_scenario = D[3:20,,], y = y)
tree = build_tree(data = flare_data, train_ind = train_samples)
print(tree,"class0_count","class1_count")
plot(tree)


# test the tree
flare_test = function(data,tree,test_ind){
  # data ratio
  train_ind = setdiff(1:length(data$y),test_ind)
  y_train = data$y[train_ind]
  #ratio = length(which(y_train==1))/length(y_train)
  
  # if it reaches a leaf node, stop and output the majority class at the leaf node
  #if (tree$leafCount==1){return(ifelse(tree$class1_count/(tree$class0_count+tree$class1_count)>ratio,1,0))} 
  
  if (tree$leafCount==1){return(tree$majority)}
  
  # otherwise we have to decide where to go
  w = tree$weight
  medoid_0 = as.integer(tree$children[[1]]$name)
  medoid_1 = as.integer(tree$children[[2]]$name)
  #p = ncol(data) - 1
  #x0 = data[medoid_0,1:p]
  #x1 = data[medoid_1,1:p]
  #feature = matrix(feature,ncol = 1)
  
  dist_mat = dist_cal(weight = w, dist_data = data$dist_scenario)
  d0 = dist_mat[medoid_0,test_ind]
  d1 = dist_mat[medoid_1,test_ind]
  
  if (d0<d1){
    return(flare_test(data = data, tree = tree$children[[1]], test_ind))
  }else{
    return(flare_test(data = data, tree = tree$children[[2]], test_ind))
  }
}

test_index = matrix(setdiff(1:length(y),train_samples),ncol = 1)
#test_index = matrix(train_samples,ncol = 1)
flare_pred = apply(test_index,1,function(x) flare_test(data = flare_data, tree = tree, test_ind = x))
acc = sum(flare_pred == y[test_index])/length(flare_pred)
acc


# visualize some of the time series
traindata = np$load("flare_train.npy")
feature_dim = which(tree$weight==max(tree$weight))+2

f1 = traindata[153,,feature_dim]
f2 = traindata[12,,feature_dim]
f3 = traindata[192,,feature_dim]
f4 = traindata[88,,feature_dim]

plotdata = data.frame(f1=f1,f2=f2,f3=f3,f4=f4,t=1:30)

library(ggplot2)
gp = ggplot(data = plotdata)+
  geom_line(aes(x = t, y = f1,color = "153"))+
  geom_line(aes(x = t, y = f2,color = "12"))+
  geom_line(aes(x = t, y = f3,color = "192"))+
  geom_line(aes(x = t, y = f4,color = "88"))

gp




# test on the log transformed data

D_log = np$load("Spar_dist_log.npy")
flare_data_log = list(dist_scenario = D_log[3:20,,], y = y)
tree_log = build_tree(data = flare_data_log, train_ind = train_samples)
plot(tree_log)
print(tree_log,"class0_count","class1_count")

test_index = matrix(setdiff(1:length(y),train_samples),ncol = 1)
#test_index = matrix(train_samples,ncol = 1)
flare_pred = apply(test_index,1,function(x) flare_test(data = flare_data_log, tree = tree_log, test_ind = x))
acc = sum(flare_pred == y[test_index])/length(flare_pred)
acc



# visualize some features
traindata = np$load("flare_train_log.npy")
node = tree_log
# +2 because the time series data contain two extra features in front of all other features, and the two features are not useful for flare prediction
feature_dim = c(5,6,11,12)+2
num_feature = length(feature_dim)
s = c(151,171,32,59)

fdata = NULL
for (i in 1:num_feature){
  for (j in s){
    fdata = c(fdata,traindata[j,,feature_dim[i]])
  }
}

t = 1:30

sclass = paste(" (Class ",as.character(y_train[s]),")",sep = "")
sclass = paste0(s,sclass,seq = "")
sclass = rep(sclass,num_feature)
fname = feature_name[feature_dim-2]

plotdata = data.frame(feature = fdata, sample = as.factor(rep(sclass,each = length(t))), 
                      time = rep(t, 4*num_feature), label = rep(y_train[s],each = length(t)),
                      feature_name = rep(fname,each = 4*length(t)))


library(ggplot2)
gp = ggplot(data = plotdata,group = feature_name)+
  geom_line(aes(x = time, y = feature, color = sample))+
  facet_wrap(~feature_name)+
  theme_bw()
  

gp






# # plot the points based on the distance matrix using multidimensional scaling
# dist_mat = dist_cal(weight = tree$weight, dist_data = D[3:20,,])
# dist_train = dist_mat
# fit = cmdscale(d = dist_train, k = 2, eig = T)
# x = fit$points[,1]
# y = fit$points[,2]
# label = c(y_train,y_test)
# 
# multD = data.frame(x = x, y = y ,label = as.factor(label))
# 
# Dp = ggplot(data = multD, aes(x = x, y = y, group = label))+
#   geom_point(aes(color = label))
# 
# Dp


# visualize the tree in network
library(networkD3)
tree_log_net = ToListExplicit(x = tree_log, unname = T)
radialNetwork(tree_log_net)

# visualize feature weight
feature_weight = data.frame(weight = c(tree_log$weight, tree_log$`171`$weight, tree_log$`151`$weight),
                            feature_name = rep(feature_name,3), node = rep(c("root","node171","node 151"),
                                                                           each = length(feature_name)))

weight_plot = ggplot(data = feature_weight,aes(x = feature_name, y = weight, fill = node))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  ylab("Feature Weight")+
  xlab("Feature Name")+
  theme(axis.text.x = element_text(angle = 90))

weight_plot








