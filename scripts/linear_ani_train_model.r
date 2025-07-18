setwd("/mnt/sdb/data/kssd3ani_project/ANIu_vs_kssd3/same_species")
#read in data
#dataset1<-read.table("Neisseria_meningitidis14.pkssd14_vs_ANIu")
#dataset2<-read.table("Kingella_kingae14.pkssd14_vs_ANIu")
Nayfach<-read.table("Nayfach52k14.pkssd14_vs_ANIm")
nrow(Nayfach)
Nayfach<-Nayfach[,c(3:10)]
header=c("X_size","N_conflict", "Y_size","XnY_ctx","N_diff_obj","N_diff_obj_section","N_mut2_ctx","ANIm")
names(Nayfach)=header
Nayfach$y=1-Nayfach$ANIm
#Nayfach$y= ifelse(Nayfach$ANIm <=1 , 1-Nayfach$ANIm, 1-Nayfach$ANIm/100)  
# Define the objective function to minimize (e.g., residual standard error)
simple_objective_fn <- function(par, data) {
  denom <- 1 / (par[1] * data$XnY_ctx + par[2] * data$N_diff_obj_section + par[3] * data$N_mut2_ctx + par[4]*data$N_diff_obj)
  model <- lm(y ~ (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom, data = transform(data, denominator = denom))
  return(summary(model)$sigma)  # Residual standard error
}
# Initial guess for a, b, c
simple_start_par <- c(1, 1, 1, 1)
# Run optimization
simple_opt_result <- optim(par = simple_start_par, fn = simple_objective_fn, data = Nayfach)
# Best coefficients
simple_opt_result$par
Nayfach$denominator <- 1 / (simple_opt_result$par[1]* Nayfach$XnY_ctx + simple_opt_result$par[2] * Nayfach$N_diff_obj_section 
                            + simple_opt_result$par[3] * Nayfach$N_mut2_ctx + simple_opt_result$par[4]*Nayfach$N_diff_obj )

simple_lm <- lm(y ~ ( XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denominator, data = Nayfach)
summary(simple_lm)
cor(predict(simple_lm),Nayfach$y)
plot(predict(simple_lm),Nayfach$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

#multi
# -------- Objective Function with 3 denominators, each using 4 terms --------
objective_fn <- function(par, data) {
  # Unpack parameters
  a1 = par[1]; b1 = par[2]; c1 = par[3]; d1 = par[4]
  a2 = par[5]; b2 = par[6]; c2 = par[7]; d2 = par[8]
  a3 = par[9]; b3 = par[10]; c3 = par[11]; d3 = par[12]
  epsilon <- 1e-8  # Small number to avoid division by 0
  # Compute denominators
  data$denom1 <- 1 / (a1 * data$XnY_ctx + b1 * data$N_diff_obj_section + c1 * data$N_mut2_ctx + d1 * data$N_diff_obj + epsilon)
  data$denom2 <- 1 / (a2 * data$XnY_ctx + b2 * data$N_diff_obj_section + c2 * data$N_mut2_ctx + d2 * data$N_diff_obj + epsilon)
  data$denom3 <- 1 / (a3 * data$XnY_ctx + b3 * data$N_diff_obj_section + c3 * data$N_mut2_ctx + d3 * data$N_diff_obj + epsilon)
  # Fit linear model with all denominators
  fit <- lm(y ~ (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom1 +
              (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom2 +
              (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom3,
            data = data)
  # Objective: residual standard error
  return(summary(fit)$sigma)
}
# Starting values for 12 parameters (a1,b1,c1,d1,...a3,b3,c3,d3)
init_params <- rep(1, 12)
# Optimize
opt_result <- optim(
  par = init_params,
  fn = objective_fn,
  data = Nayfach,
  method = "Nelder-Mead",
  control = list(maxit = 1000, reltol = 1e-6)
)

# View optimal parameters
cat("Optimal parameters:\n")
param_names <- paste0(rep(c("a","b","c","d"), 3), rep(1:3, each = 4))
print(setNames(opt_result$par, param_names))
# Extract optimized parameters
p <- opt_result$par
# Compute final denominators
epsilon <- 1e-8
Nayfach$denom1 <- 1 / (p[1] * Nayfach$XnY_ctx + p[2] * Nayfach$N_diff_obj_section + p[3] * Nayfach$N_mut2_ctx + p[4] * Nayfach$N_diff_obj + epsilon)
Nayfach$denom2 <- 1 / (p[5] * Nayfach$XnY_ctx + p[6] * Nayfach$N_diff_obj_section + p[7] * Nayfach$N_mut2_ctx + p[8] * Nayfach$N_diff_obj + epsilon)
Nayfach$denom3 <- 1 / (p[9] * Nayfach$XnY_ctx + p[10] * Nayfach$N_diff_obj_section + p[11] * Nayfach$N_mut2_ctx + p[12] * Nayfach$N_diff_obj + epsilon)
# Final model using optimized denominators
final_model <- lm(y ~ (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom1 +
                    (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom2 +
                    (XnY_ctx + N_diff_obj_section + N_mut2_ctx + N_diff_obj) * denom3,
                  data = Nayfach)
# Show summary
summary(final_model)

cor(predict(final_model),Nayfach$y)
plot(predict(final_model),Nayfach$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

#####other features(may not be used)

Nayfach$Min_size=ifelse(Nayfach$X_size < Nayfach$Y_size,Nayfach$X_size,Nayfach$Y_size)
Nayfach$Jcd = Nayfach$XnY_ctx/(Nayfach$X_size + Nayfach$Y_size - Nayfach$XnY_ctx)
Nayfach$Md = -log(2*Nayfach$Jcd/(1+Nayfach$Jcd))
Nayfach$Aaf= -log(Nayfach$XnY_ctx/Nayfach$Min_size)
Nayfach$cod=Nayfach$N_diff_obj/Nayfach$XnY_ctx

Nayfach$r1= ifelse(Nayfach$y < 0.01,1,Nayfach$y*9*Nayfach$XnY_ctx/Nayfach$N_diff_obj_section)

Nayfach$ratio=(Nayfach$N_mut2_ctx+1)/(Nayfach$N_diff_obj-Nayfach$N_mut2_ctx+1)
Nayfach$ratio2=(Nayfach$N_diff_obj_section-Nayfach$N_diff_obj+1)/(Nayfach$N_diff_obj+1)

Nayfach$y2x=ifelse(Nayfach$cod < 0.01,1,Nayfach$y/Nayfach$cod)
Nayfach$r2c=(Nayfach$ratio+0.001)/(Nayfach$cod+0.001)
Nayfach$r2c2=(Nayfach$ratio2+0.001)/(Nayfach$cod+0.001)
Nayfach$diff1 = Nayfach$N_diff_obj-Nayfach$N_mut2_ctx
Nayfach$diff2 = Nayfach$N_diff_obj_section-Nayfach$diff1
#Nayfach$r3= (Nayfach$diff1+1)/(Nayfach$diff2+1)
Nayfach$r3= (Nayfach$diff1+1)/(Nayfach$N_diff_obj_section+1)
Nayfach$best<-(Nayfach$N_diff_obj_section+Nayfach$N_mut2_ctx)*Nayfach$denominator

#complicate linear model
m2lm<-lm(y~  best*XnY_ctx*N_mut2_ctx*N_diff_obj*r3, data=Nayfach)
summary(m2lm)
cor(predict(m2lm),Nayfach$y)
plot(predict(m2lm),Nayfach$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
###test_data
test<-read.table("Kingella_kingae14.pkssd14_vs_ANIu")
test<-read.table("Neisseria_meningitidis14.pkssd14_vs_ANIu")
test<-test[,c(3:10)]
names(test)=header

test$denominator <- 1 / (simple_opt_result$par[1]* test$XnY_ctx + simple_opt_result$par[2] * test$N_diff_obj_section 
                            + simple_opt_result$par[3] * test$N_mut2_ctx + simple_opt_result$par[4]*test$N_diff_obj )
test$y=1-test$ANIm/100  


test$denom1 <- 1 / (p[1] * test$XnY_ctx + p[2] * test$N_diff_obj_section + p[3] * test$N_mut2_ctx + p[4] * test$N_diff_obj + epsilon)
test$denom2 <- 1 / (p[5] * test$XnY_ctx + p[6] * test$N_diff_obj_section + p[7] * test$N_mut2_ctx + p[8] * test$N_diff_obj + epsilon)
test$denom3 <- 1 / (p[9] * test$XnY_ctx + p[10] * test$N_diff_obj_section + p[11] * test$N_mut2_ctx + p[12] * test$N_diff_obj + epsilon)

final_predictions <- predict(final_model, newdata = test)

cor(final_predictions,test$y)
plot(final_predictions,test$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

simple_predictions <- predict(simple_lm, newdata = test)

cor(simple_predictions,test$y)
plot(simple_predictions,test$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

test$diff1 = test$N_diff_obj-test$N_mut2_ctx
test$diff2 = test$N_diff_obj_section-test$diff1
test$r3= (test$diff1+1)/(test$N_diff_obj_section+1)

test$Min_size=ifelse(test$X_size < test$Y_size,test$X_size,test$Y_size)
test$ratio=(test$N_mut2_ctx+1)/(test$N_diff_obj-test$N_mut2_ctx+1)
test$ratio2=(test$N_diff_obj_section-test$N_diff_obj+1)/(test$N_diff_obj+1)

test$Jcd= test$XnY_ctx/(test$X_size + test$Y_size - test$XnY_ctx)
test$Md = -log(2*test$Jcd/(1+test$Jcd))
test$Aaf= -log(test$XnY_ctx/test$Min_size)
test$cod=test$N_diff_obj/test$XnY_ctx
test$best<-(test$N_diff_obj_section+test$N_mut2_ctx)*test$denominator
predictions <- predict(m2lm, newdata = test)

cor(predictions,test$y)

plot(predictions,test$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

#########old code
##
Cl<-read.table("Clostridioides_difficile5.pkssd3_vs_ANIu")
train.Cl<- Cl[sample(nrow(Cl), 500), ]
Nm<-read.table("Neisseria_meningitidis5.pkssd3_vs_ANIu")
train.Nm <- Nm[sample(nrow(Nm), 5000), ]
Ecoli<-read.table("Escherichia_coli5.pkssd3_vs_ANIu")
train.Ecoli<-Ecoli[sample(nrow(Ecoli), 500), ]
Kig<-read.table("Kingella_kingae5.pkssd3_vs_ANIu")
train.Kig<-Kig[sample(nrow(Kig), 500), ]

Vpara<-read.table("Vibrio_parahaemolyticus5.pkssd3_vs_ANIu")
train.Vpara<-Vpara[sample(nrow(Vpara), 500), ]
Cbuty<-read.table("Clostridium_butyricum5.pkssd3_vs_ANIu")
train.Cbuty<-Cbuty[sample(nrow(Cbuty),500),]
Bmel<-read.table("Brucella_melitensis5.pkssd3_vs_ANIu")
train.Bmel<-Bmel[sample(nrow(Bmel),500),]
PEpara<-read.table("Pseudomonas_E_paracarnis5.pkssd3_vs_ANIu")
train.PEpara<-PEpara[sample(nrow(PEpara),500),]

Nayfach<-read.table("Nayfach52k.pkssd3_vs_ANIm")
Nayfach[,ncol(Nayfach)]=Nayfach[,ncol(Nayfach)]*100
train.Nayfach<-Nayfach[sample(nrow(Nayfach),1000),]
train.Nayfach<-train.Nayfach[train.Nayfach[,ncol(Nayfach)] > 95,]
#Nayfach11<-read.table("Nayfach52k11.pkssd11_vs_ANIm")
#Nayfach11[,ncol(Nayfach11)]= Nayfach11[,ncol(Nayfach11)]*100
#Nayfach<-Nayfach[sample(nrow(Nayfach),1000),]
#lt95Nayfach<-read.table("lt095.Nayfach52k5.pkssd3_vs_ANIm")
#lt95Nayfach[,ncol(lt95Nayfach)]=lt95Nayfach[,ncol(lt95Nayfach)]*100
#lt95Nayfach<-lt95Nayfach[sample(nrow(lt95Nayfach),1000),]

merge<-rbind(Nayfach)
#merge<-rbind(train.Nayfach,train.Cl,train.Nm,train.Ecoli,train.Kig,
#             train.Vpara,train.Cbuty,train.Bmel,train.PEpara,train.Kleb)
subcols<-c(3,5,6,7,8,13)
subcols.name=c("X_size","Y_size","XnY_ctx","N_diff_obj","N_diff_obj_section","ANIu")
merge<-merge[,subcols]
colnames(merge)<-subcols.name

merge$y=1-merge$ANIu/100
#merge$Min_size= ifelse(merge$X_size<merge$Y_size,merge$X_size,merge$Y_size)
#mlm<-lm(y~XnY_ctx*N_diff_obj*N_diff_obj_section*I(1/XnY_ctx)*I(1/Min_size),data=merge)
mlm<-lm(y~XnY_ctx*N_diff_obj*N_diff_obj_section*I(1/XnY_ctx),data=merge)
summary(mlm)
cor(predict(mlm),merge$y)
plot(predict(mlm),merge$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

#test<-read.table("Neisseria_meningitidis5.pkssd3_vs_ANIu")
test<-rbind(Cl,Nm,Ecoli,Kig,Vpara,Cbuty,Bmel,PEpara,Kleb)
test<-Kig
test<-test[,subcols]
colnames(test)<-subcols.name
test$Min_size= ifelse(test$X_size<test$Y_size,test$X_size,test$Y_size)
test$y=1-test$ANIu/100
#
predictions <- predict(mlm, newdata = test)
cor(test$N_diff_obj/test$XnY_ctx,test$y)
cor(predictions,test$y)
plot(predictions,test$y)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)


###
Nayfach13<-read.table("Nayfach52k13.pkssd13_vs_ANIm")
Nayfach13<-Nayfach13[sample(nrow(Nayfach13),10000),]

Nayfach13<-Nayfach13[,c(3:10)]
names(Nayfach13)=c("X_size","N_conflict", "Y_size","XnY18","N_diff_obj18","XnY28","N_diff_obj28","ANIm")
Nayfach13$y=   1-Nayfach13$ANIm
Nayfach13$dist1 = Nayfach13$N_diff_obj18/Nayfach13$XnY18
Nayfach13$dist2= Nayfach13$N_diff_obj28/Nayfach13$XnY28
Nayfach13$rate= (Nayfach13$dist2+0.0001)/(Nayfach13$dist1+0.0001)
Nayfach13$Min_size = ifelse(Nayfach13$X_size < Nayfach13$Y_size,Nayfach13$X_size,Nayfach13$Y_size)
Nlm<-lm(y~N_diff_obj18*N_diff_obj28*XnY18*I(1/XnY18)*I(1/Min_size) ,data=Nayfach13)
summary(Nlm)
cor(predict(Nlm),Nayfach13$y)
plot(predict(Nlm),Nayfach13$y)

test<-read.table("Neisseria_meningitidis14.pkssd14_vs_ANIu")
test2<-read.table("Kingella_kingae14.pkssd14_vs_ANIu")
test$V10<-test$V10/100
test2$V10<-test2$V10/100

