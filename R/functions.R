#' @title Convolution function
#' @description Simple convolution function is to extract features with using filters(kernel),here only consider the case where the number of channels is 1, that is, to process gray images
#' @param img the matrix of picture pixel
#' @param kernel_size is given according to the experimental requirements,and kernel, the filter is generated randomly
#' @return the matrix after convolution,the size is ((img_nrow - kernel_row)/stride + 1)*((img_col - kernel_col)/stride + 1)
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' simple_conv2d(img = matrix(data = rnorm(100), ncol = 10), kernel_size = c(4, 4))
#' }
#' @export
simple_conv2d <- function(img, kernel_size) {
  img_col = ncol(img)
  img_nrow = nrow(img)
  kernel_col = kernel_size[1]
  kernel_row = kernel_size[2]
  #Let stride = 1
  #kernel is generated randomly 
  kernel <- matrix(data = rnorm(n = kernel_col * kernel_row, mean = 0, sd = 1), ncol = kernel_col)
  
  result <- matrix(nrow = (img_nrow - kernel_row + 1), ncol = (img_col - kernel_col + 1))
  
  for (i in c(1:(img_col - kernel_col+1))) {
    for (j in c(1:(img_nrow - kernel_row +1))) {
      result[j, i] <- sum(img[j:(j+kernel_row-1), i:(i+kernel_col-1)] * kernel)
    }
  }
  
  return(result)
}


#' @title Max Pooling function
#' @description This function is to reduce dimensionality and prevent overfitting
#' @param hidden the matrix of the neurons after the activation function
#' @return the matrix after Pooling,the size is half of the previous
#' @examples
#' \dontrun{
#' MaxPool2d(hidden = matrix(data = rnorm(100), ncol = 10))
#' }
#' @export
MaxPool2d <- function(hidden) {
  hidden_col = ncol(hidden)
  hidden_nrow = nrow(hidden)
  #the kernel with a size of 2*2 and a stride of 1 is usually selected
  result <- matrix(nrow = ((hidden_nrow - 2)/2 + 1), ncol = ((hidden_col - 2)/2 + 1))
  
  for (i in c(1:((hidden_col - 2)/2 + 1))) {
    for (j in c(1:((hidden_nrow - 2)/2 + 1))) {
      result[j, i] <- max(hidden[j:((j+2)/2-1), i:(i+2-1)] )
    }
  }
  return(result)
}

