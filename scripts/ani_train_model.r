#library(gbm)
library(xgboost)
library(h2o)
setwd("/mnt/sdb/data/kssd3ani_project/ANIm")
sk<-read.table("dedup.out_Nayfach52k_kssd3f8C9O7_vs_ANIm_skani.out",header=T)
cor(sk$ANIm,sk$Skani_ANI)
#ANIm<-sk[,c(4,5,6,7,8,13)]

ANIm<-read.table("cut-f1-8_14.dedup.out_Nayfach52k_kssd3f8C9O7_vs_ANIm.out",header=T)
#ANIm<-ANIm[ANIm$ANI > 98 & ANIm$Qry_align_fraction >0.35 & ANIm$Ref_align_fraction >0.35,]
ANIm$ANIm=ANIm$ANIm*100
ANIm <- ANIm[, -c(1,2,8)]
ANIm<-ANIm[sample(1:nrow(ANIm),  90000),]

#ANIm<-ANIm[ANIm$Qry_align_fraction > 0.3 & ANIm$Ref_align_fraction >0.3 & ANIm$learned_ANI >96 ,]
#plot(ANIm$ANI,ANIm$ANIm)
train_idx <- sample(1:nrow(ANIm),  0.8 * nrow(ANIm))
train_data <- ANIm[train_idx, ]
test_data <- ANIm[-train_idx, ]

set.seed(123)
split_idx <- sample(nrow(train_data), 0.8 * nrow(train_data))  # 80%训练，20%验证$
train_sub <- train_data[split_idx, ]
valid_sub <- train_data[-split_idx, ]

# 转换为 DMatrix 格式
dtrain <- xgb.DMatrix(data = as.matrix(train_sub[, -ncol(train_sub)]), label = train_sub$ANIm)
dvalid <- xgb.DMatrix(data = as.matrix(valid_sub[, -ncol(valid_sub)]), label = valid_sub$ANIm)
feature_names <- colnames(dtrain)  # replace with your actual DMatrix variable
# 定义 watchlist
watchlist <- list(train = dtrain, eval = dvalid) 
xgb_model <- xgb.train(
  data = dtrain,
  objective = "reg:squarederror",
  max_depth = 4,
  eta = 0.1,
  nrounds = 200,
  watchlist = watchlist,          # 关键：传递验证集
  early_stopping_rounds = 10,     # 验证集性能连续10轮未提升则停止
  eval_metric = "rmse"            # 回归任务常用指标（可选：mae, mape等）
)
importance_matrix <- xgb.importance(feature_names = feature_names, model = xgb_model)
#xgb.plot.importance(importance_matrix)
# Predictions
#gbm_pred <- predict(gbm_model, test_data, n.trees = best_iter)
#cor(gbm_pred,test_data$ANIm)
cor(test_data$ANI,test_data$ANIm)
#test1<-read.table("/mnt/sdb/data/kssd3ani_project/ANIu_vs_kssd3/same_species/merge5species.kssd3f8C9O7_f8C11O5_vs_ANIu_skani")
#test1<-test1[,c(1:5,8)]
#colnames(test1)<-c(feature_names,"ANIu")
#xgb_pred <- predict(xgb_model, xgb.DMatrix(data = as.matrix(test1[, -ncol(test1)])))
#cor(xgb_pred,test1$ANIu)
loaded_model <- xgb.load("/home/ubt/work1/tools/KSSD3/ani_models/f8C9O7_ANI98.xgb")
xgb_pred <- predict(xgb_model, xgb.DMatrix(data = as.matrix(test_data[, -ncol(test_data)])))
cor(xgb_pred,test_data$ANIm)
#tmp<-cbind(xgb_pred,test_data)
#tmp<-tmp[tmp$Qry_align_fraction>0.3 & tmp$Ref_align_fraction>0.3,]
#plot(tmp$co.distance,tmp$ANIm)
xgb.save(xgb_model, "model2.xgb")
#xgb_pred <- predict(xgb_model, xgb.DMatrix(test_data[, -ncol(test_data)]))
#summary(gbm_model, plotit = TRUE)  
# Calculate RMSE/MAE
#rmse <- sqrt(mean((test_data$Sale_Price - gbm_pred)^2))
mae <- mean(abs(test_data$Sale_Price - xgb_pred))

#load
model <- xgb.load("f8C9O7_model.xgb")

# Dump model as text
dump <- xgb.dump(model, with_stats = TRUE)
cat(dump[1:10], sep = "\n")  # Show first few lines

# Extract lines that contain features (e.g., [f4<...])
feature_lines <- grep("\\[f[0-9]+<", dump, value = TRUE)

# Extract feature indices using regex
feature_ids <- as.integer(gsub(".*\\[f([0-9]+)<.*", "\\1", feature_lines))

# Get unique features and compute max index + 1
num_features <- max(feature_ids) + 1

cat("Number of features used by the model:", num_features, "\n")



library(gbm)

# 定义超参数网格
hyper_grid <- expand.grid(
  n.trees = c(500, 1000),
  interaction.depth = c(3, 5),
  shrinkage = c(0.01, 0.1)
)

best_mse <- Inf
best_params <- list()

for (i in 1:nrow(hyper_grid)) {
  params <- hyper_grid[i, ]
  
  # 核心模型调用
  model <-  gbm(
    V6 ~ ., 
    data = train_data,
    distribution = "gaussian",  # For regression
    n.trees = 500,             # Number of trees
    interaction.depth = 4,     # Tree depth
    shrinkage = 0.01,          # Learning rate
    cv.folds = 5               # Cross-validation
  )
  
  # 获取最优树数量对应的MSE
  best_iter <- gbm.perf(model, method = "cv")
  mse <- model$cv.error[best_iter]
  
  # 更新最佳参数
  if (mse < best_mse) {
    best_mse <- mse
    best_params <- params
  }
}

# 输出最佳参数组合
print(paste("Best MSE:", round(best_mse, 4)))
print("Optimal parameters:")
print(best_params)

# Predictions
gbm_pred <- predict(gbm_model, test_data, n.trees = best_iter)
xgb_pred <- predict(xgb_model, xgb.DMatrix(test_data[, -ncol(test_data)]))

# Calculate RMSE/MAE
rmse <- sqrt(mean((test_data$Sale_Price - gbm_pred)^2))
mae <- mean(abs(test_data$Sale_Price - xgb_pred))




plot(Clod[,6],Clod[,14])



setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANI94_98/try2_gtdb_genomes_reps_r220_f8O5")
file<-read.table("ACGTonly.even_distribut1k.ANI96_98.sorted.try2_gtdb_genomes_reps_r220_f8C9O7_f8C11O5.ANIu")
plot(file[,5],file[,9])
m<-lm(V9~V1+V2+V3+V4+V5,data=file)
summary(m)
predictions <- predict(m, newdata = file)
plot(file[,5],file[,9])
plot(predictions,file[,9])
file<-file[file[,9] > 0,]
file2<-read.table("tmp.shuf")
plot(file[,6],file[,9])
train<-read.table("h0.4.tmp.shuf")
test<-read.table("t0.4.tmp.shuf")
#plot(train[,6],train[,9])
m<-lm(V9~V1+V2+V3+V4+V5,data=train)
predictions <- predict(m, newdata = test)
cor(test[,9],predictions)
plot(test[,9],predictions)

summary(m)


setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANI94_98/Klebsiella_pneumoniae2quasipneumoniae")
tmp<-read.table("noambigous.out1k.Klebsiella_pneumoniae2quasipneumoniae.f8C9O7_f8C11O5_vs_ANIu")
newd<-tmp[,c(1:6,12:19)]
predictions <- predict(m, newdata = newd)
summary(lm(newd[,9] ~ newd[,5]))
colnames(newd)<-head(colnames(tmp),14)
plot(tmp[,6],tmp[,14])
setwd("/mnt/sdb/data/ANIu_vs_kssd3/same_species/")
Clod<-read.table("53noAmbigu.Clostridioides_difficile.kssd3f8C9O7_f8C11O5_vs_ANIu_skani")
plot(Clod[,5],Clod[,13])
clod<-read.table("tmp")
new<-read.table("noambiguosnt.subset100_same_species.kssd3f8C9O7_f8C11O5_vs_ANIu")
l<-lm( V13 ~ V5 + V3 + V4 + V1 + V2 , data=tmp)
summary(l)
predictions <- predict(l, newdata = new)
cor(predictions,new[,13])
plot(predictions,new[,13])
plot(new[,5],new[,13])
summary(l)
cor(clod[,13],predict(l,val[,5]))
plot(clod[,13],predict(l))
cor(clod[,13],clod[,19])
summary(lm(clod[,13] ~ clod[,19] ))
p<-100-(100-clod[,5])*1.52
#clod<-clod[clod[,5] > 98,]
#p<-(1-clod[,4])**(1/10)
#plot(p,clod[,8])
#cor(p,clod[,8])
cor(clod[,19],clod[,13])
cor(clod[,14],clod[,8])
plot(clod[,10],clod[,13])
plot(clod[,14],clod[,8])
summary(lm(clod[,8] ~ clod[,5]))
summary(lm(clod[,8] ~ clod[,14]))
plot(p,clod[,8])
summary(lm(clod[,8] ~ p))
sum((clod[,8] - p)**2)
sum((clod[,8] - clod[,14])**2)

setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANI94_98/Klebsiella_pneumoniae2quasipneumoniae")
subk<-read.table("tmp")
plot(subk[,4],subk[,5])
plot(100-subk[,5],100-subk[,8])
a<-100-subk[,5]
b<-100-subk[,8]
cor(a,b)
summary(lm(b ~  a))
summary(lm(100 - subk[,8] ~ 100 - subk[,5]))
plot(1-subk[,4],(subk[,8]/100)**10)

skani<-read.table("subset100_same_species.skani_vs_oau")
cor(skani[,3],skani[,8])
pkssd<-read.table("subset100_same_species.perlkssd3_vs_oau")

plot(pkssd[,10],pkssd[,15])
cor(pkssd[,10],pkssd[,15])

setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANI94_98/try2_gtdb_genomes_reps_r220_f8O5")
f7C9O1<-read.table("kssd3outf7C9O1_10k.ANI94_98_vs_ANIu")
plot(f7C9O1[,5],f7C9O1[,8])
cor(f7C9O1[,5],f7C9O1[,8])

f7C9O7<-read.table("kssd3outf7C9O710k.ANI94_98_vs_ANIu")
f7C9O7<-f7C9O7[f7C9O7[,2]>0.2 & f7C9O7[,3] >0.2,]
plot(f7C9O7[,5],f7C9O7[,8])
cor(f7C9O7[,5],f7C9O7[,8])
summary(lm(f7C9O7[,8]~f7C9O7[,5]+f7C9O7[,2]+f7C9O7[,4]))

f7C11O5<-read.table("kssd3outf7C9O510k.ANI94_98_vs_ANIu")
plot(f7C11O5[,5],f7C11O5[,8])
cor(f7C11O5[,5],f7C11O5[,8])
sk<-read.table("10kskani_vs_ANIu")
plot(sk[,1],sk[,2])
cor(sk[,1],sk[,2])
j<-read.table("kssd3_vs_ANIu.10keven")
j<-j[j[,8]>60 & j[,2] > 0.2 & j[,3] > 0.2 & j[,1] > 1000,]
plot(j[,5],j[,8])
cor(j[,5],j[,8])
summary(lm(j[,8] ~ j[,5]))

j9698<-read.table("ACGTonly.even_distribut1k.ANI96_98.sorted.try2_gtdb_genomes_reps_r220.f8C11O5_vs_ANIu")
j9698<-j9698[j9698[,8]>88,]
plot(j9698[,5],j9698[,8])
cor(j9698[,5],j9698[,8])
summary(lm(j9698[,8] ~ j9698[,5]))
setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANI94_98/Klebsiella_pneumoniae2quasipneumoniae")
kleb<-read.table("out1k.Klebsiella_pneumoniae2quasipneumoniae.f7C11O5_vs_ANIu")
kleb<-kleb[kleb[,8]>93.6 & kleb[,8] < 94.6,]
plot(kleb[,5],kleb[,8])
cor(kleb[,5],kleb[,8])
setwd("/mnt/sdb/data/ANIu_batch2/new")
k<-read.table("test1")
k<-k[k[,1] >0.3 & k[,3]>93.5 &  k[,3] < 94.8 ,]

plot(k[,2],k[,3])
cor(k[,2],k[,3])
summary(lm(k[,3]~k[,2]))

setwd("/mnt/sdb/data/ANIu_vs_kssd3/ANIu_batch2/batch101_300")
t<-read.table("batch101_300.ANIuVSKSSD3f8C11O5")

u<-read.table("out1k.sample_even_distribution")
hist(u[,7],breaks=100)
t<-t[t[,3] >90,]
plot(t[,3],t[,4])

kd1<-c(subk[,5],t[,4])
kd2<-c(subk[,8],t[,3])
plot(kd2,kd1)
summary(lm(kd2~kd1))

aniu<-c(t[,3],k[,3])
anik<-c(t[,4],k[,2])
cor(aniu,anik)
plot(aniu,anik)
m<-t
min(t[,4])
m[,3] = 100-m[,3]
m[,4] = 100-m[,4]
plot(m[,3],m[,4])
plot(t[,3],t[,4])

summary(lm(m[,3]~m[,4]))
summary(lm(t[,3]~t[,4]))

cor(m[,3],m[,4])

t8<-read.table("batch101_300.ANIuVSKSSD3f8C8O5")
t8<-t8[t8[,4] > 0,]
t8<-t8[t8[,3] > 95,]
plot(t8[,3],t8[,4])
summary(lm(t8[,3]~t8[,4]))

t9<-read.table("batch101_300.ANIuVSKSSD3f8C9O5")
t9<-t9[t9[,4] > 0,]
t9<-t9[t9[,3] > 95,]
summary(lm(t9[,3]~t9[,4]))

plot(t9[,3],t9[,4])
cor(t9[,3],t9[,4])

t10<-read.table("batch101_300.ANIuVSKSSD3f8C10O5")
t10<-t10[t10[,4] > 0,]
t10<-t10[t10[,3] > 95,]
summary(lm(t10[,3]~t10[,4]))
plot(t10[,3],t10[,4])
cor(t10[,3],t10[,4])
