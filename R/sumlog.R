#'@title  Internal function for EM
#'@description  Internal function for EM
#'@param A1 an integer
#'@param A2 an integer


sumlog <- function (A1, A2)
{

  if(A1 > A2){
    res = A1 + log(1 + exp(A2 - A1))
  }else{
    res = A2 + log(exp(A1 - A2) + 1)
  }

  return (res)

}
