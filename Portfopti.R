## install.packages("Rglpk")

library(Rglpk)
Return_data <- readr::read_rds("E:/RWD/CDaRPortf/return_data.rds")
alpha <- 0.8
omega <- 0.2
Portfopti <- function(Return_data,alpha,omega){
  ## Inputs: Returndata---The logarithm return data
  ## alpha----The confident Level
  ## omega----The up bound of CDaR
  
data_dim <- dim(Return_data)
sum_data <- Return_data
 for(i in 2:data_dim[1]){
   sum_data[i,] <- apply(Return_data[1:i,],2,sum)
 }
 
A <- matrix(0,data_dim[1],data_dim[1])
for(i in 1:data_dim[1]){
  A[i,i] <- -1
  A[i,i-1] <- 1
}

Opti_obj <- cbind(0,as.matrix(sum_data[data_dim[1],]), matrix(0,1,data_dim[1]),
                  matrix(0,1,data_dim[1]))
c1 <- cbind(1, matrix(0,1,data_dim[2]),
         1/((1-alpha)*data_dim[1])*matrix(1,1,data_dim[1]),matrix(0,1,data_dim[1]));
c2 <- cbind(-matrix(1,data_dim[1],1),-as.matrix(sum_data),-diag(rep(1,data_dim[1])),
                                           diag(rep(1,data_dim[1])));
c3 <- cbind(matrix(0,data_dim[1],1),as.matrix(sum_data),matrix(0,data_dim[1],data_dim[1]),
        -diag(rep(1,data_dim[1])))

c4 <- cbind(matrix(0,data_dim[1],1),matrix(0,data_dim[1],data_dim[2]),
        matrix(0,data_dim[1],data_dim[1]),A)

mat <- rbind(c1,c2,c3,c4)

dir <- matrix("<=",3*data_dim[1]+1,1)

rhs <- cbind(omega,matrix(0,1,3*data_dim[1]))
bounds <- list(lower = list(ind = c(1L),val = c(-Inf)),
              upper = list(ind = c(2:(data_dim[2]+1)),val = rep(1,data_dim[2])))

Retur <- Rglpk::Rglpk_solve_LP(Opti_obj,mat,dir,rhs,bounds,max=TRUE)

Retur$solution[2:(data_dim[2]+1)]
}


Portfopti(Return_data,0.8,0.2)

