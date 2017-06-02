#
# Source the G-test
#currdir<-getwd()
#source(paste(currdir,"/G_test.R",sep=""))
#
#
# FUNCTION: Woolf-test
# The script comes from the help page for the mantelhaen.test()
# ?mantelhaen.test()
woolf.test <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  woolf <- sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2)
  df <- k-1
  p <- 1 - pchisq(woolf, df)
  dat <- c(woolf,df,p)
  names(dat) <- c("Woolf", "df", "p-value")
  dat
}

# FUNCTION: Get GLM results from array ##
#get GLM data-set from k-way table
get_glm_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array <- array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  pop <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      pop[j] <- dimnames(array)[3][[1]][k]
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep)
  mod <- anova(glm(
    cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
    family = "binomial"),test="LRT")
  return(mod)
}
# FUNCTION: convert array to data.frame
get_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array<-array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  pop <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      pop[j] <- dimnames(array)[3][[1]][j]
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"pop"=pop)
  return(d)
}
