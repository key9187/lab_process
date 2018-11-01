rm(list=ls())
## 3. kospi data set example ---------------------------------------------------------

if(!require(MASS))install.packages("MASS");library(MASS)
setwd("C:/Users/UOS/OneDrive - 서울시립대학교/key/세미나 준비/optimal portfolio")
load("example_kospi_dat.RData")

# # 일일 수익률(rate of return) = (나중주가-처음주가)/처음주가
# example_rr <- example_kospi_dat/rbind(NA, example_kospi_dat[-nrow(example_kospi_dat),])-1
# example_rr <- example_rr[-1,]

# 일일 로그수익률(log rate of return) = log(나중주가/처음주가)
example_rr <- log(example_kospi_dat)-log(rbind(NA, example_kospi_dat[-nrow(example_kospi_dat),]))
example_rr <- example_rr[-1,]

# init value
n <- dim(example_rr)[2]; t <- dim(example_rr)[1]
mu_0 <- 10   # 목표수익률 0.1%
R <- as.matrix(example_rr)   # kospi자료 6개 기업(1~8월 data)
R_i <- ginv(R)
one <- rep(mu_0,t)

rho <- 1   # 임의의 값 지정 
mu_hat <- colMeans(R)*100   # 알려진 수익률 추정해야 함

lambda <- 10   # 임의의 값 지정
A <- rbind(diag(x=n),1,t(mu_hat))
B <- rbind(-diag(x=n), 0, 0)
c1 <- c(rep(0,n), 1, mu_0)
z_init <- c(rep(0,n))
u_init <- c(rep(1,n+2))


# x update

tmp_inverse <- ginv(t(R) %*% R + rho * t(A) %*% A )
tmp_tRone <- t(R) %*% one 

func_x <- function(z, u, tmp_inverse, tmp_tRone, rho, A, B, c1){
  x <- tmp_inverse %*% (tmp_tRone - rho * t(A) %*% (B %*% z - c1 + u) )
  return(x)
}  


# z update

tmp_lr <- lambda/rho                                                                    

func_z <- function(u, x, tmp_lr, A, B, c1){
  z <- ifelse(abs(t(B) %*% (A %*% x - c1 + u)) < tmp_lr, 0,
              - t(B) %*% (A %*% x - c1 + u) + sign(t(B) %*% (A %*% x - c1 + u)) * tmp_lr)
  return(z)
}


# u update

func_u <- function(x, z, u, A, B, c1){
  u <- u + (A %*% x + B %*% z - c1)
  return(u)
}


# 반복하여 구하기

i=1

while(i <= 1000){
  
  tmp_x <- func_x(z = z_init, u = u_init, tmp_inverse=tmp_inverse, tmp_tRone=tmp_tRone, rho=rho, A=A, B=B, c1=c1)
  
  tmp_z <- func_z(u = u_init, x = tmp_x, tmp_lr=tmp_lr, A=A, B=B, c1=c1)
  
  tmp_u <- func_u(x = tmp_x, z = tmp_z, u = u_init, A=A, B=B, c1=c1)

  #stop criteria
  if(any(is.na(tmp_x))) break
  if(i %% 200 == 0)
    cat('i:::', i , '\n','x:::', tmp_x, '\n', 
        'loss:::', t(one - R %*% tmp_x) %*% (one - R %*% tmp_x) + lambda * sum(abs(tmp_x)),'\n',
        '===========================','\n')
  
  z_init <- tmp_z
  u_init <- tmp_u
  
  loss[i] <- t(one - R %*% tmp_x) %*% (one - R %*% tmp_x) + lambda * sum(abs(tmp_x))
  
  i <- i+1
}


solution <- cbind(c(tmp_x, NA, NA), c(tmp_z, NA, NA), tmp_u); colnames(solution) <- c("x", "z", "u"); solution
# 제약조건 확인
sum(tmp_x)
c(t(mu_hat) %*% tmp_x, mu_0)
# A %*% tmp_x + B %*% tmp_z - c1
# plot(A %*% tmp_x + B %*% tmp_z - c1)

# x : optimal z, u
(t(R) %*% R + rho * t(A) %*% A) %*% tmp_x - tmp_tRone + rho * t(A) %*% (B %*% tmp_z - c1 + tmp_u)

# z : optimal x, u
- rho * t(B) %*% (A %*% tmp_x - c1 + tmp_u)
c(-tmp_lr, tmp_lr)

# u : optimal z, u
rho * (A %*% tmp_x + B %*% tmp_z - c1)




