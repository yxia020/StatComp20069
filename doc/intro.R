## ----eval=FALSE---------------------------------------------------------------
#  simple_conv2d <- function(img, kernel_size) {
#    img_col = ncol(img)
#    img_nrow = nrow(img)
#    kernel_col = kernel_size[1]
#    kernel_row = kernel_size[2]
#    #generate norm
#    kernel <- matrix(data = rnorm(n = kernel_col * kernel_row, mean = 0, sd = 1), ncol = kernel_col)
#  
#    result <- matrix(nrow = (img_nrow - kernel_row + 1), ncol = (img_col - kernel_col + 1))
#  
#    for (i in c(1:(img_col - kernel_col+1))) {
#      for (j in c(1:(img_nrow - kernel_row +1))) {
#        result[j, i] <- sum(img[j:(j+kernel_row-1), i:(i+kernel_col-1)] * kernel)
#      }
#    }
#  
#    return(result)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  MaxPool2d <- function(hidden) {
#    hidden_col = ncol(hidden)
#    hidden_nrow = nrow(hidden)
#    #generate norm
#    result <- matrix(nrow = ((hidden_nrow - 2)/2 + 1), ncol = ((hidden_col - 2)/2 + 1))
#  
#    for (i in c(1:((hidden_col - 2)/2 + 1))) {
#      for (j in c(1:((hidden_nrow - 2)/2 + 1))) {
#        result[j, i] <- max(hidden[j:((j+2)/2-1), i:(i+2-1)] )
#      }
#    }
#  
#    return(result)
#  }

## ----warning=FALSE, eval=FALSE------------------------------------------------
#  res1 <- simple_conv2d(img = matrix(data = rnorm(100), ncol = 10), kernel_size = c(5, 5))
#  res2 <- MaxPool2d(hidden = matrix(data = rnorm(100), ncol = 10))
#  res1
#  res2
#  dim(res1)
#  dim(res2)

