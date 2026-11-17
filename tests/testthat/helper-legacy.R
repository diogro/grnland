# Closed-form genotype -> phenotype maps and sensitivity matrices copied
# from the original analysis scripts (GRN_Landscapes/Figs34.r, Fig2.r).
# The Fig2.r functions are named J* there but are the sensitivities dx*/dalpha.
motif_1 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- a2 * (x1 / (th + x1)) / g
  
  return(c(x1, x2))
}

motif_2 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- (a2 * (1 - ((a1 / g) / (th + ((a1) / (g)
  ))))) / (g)
  
  return(c(x1, x2))
}

motif_3 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (2 * th * a1 * sqrt(g)) / ((-a1 + th * g) * sqrt(g) + sqrt(4 * a1 *
                                                                     a2 * g + ((a1 + th * g)^2) * g))
  x2 <- (-(a1 + th * g) * sqrt(g) + sqrt(4 * a1 * a2 * g + ((a1 + th * g)^2) *
                                           g)) / (2 * g * sqrt(g))
  
  return(c(x1, x2))
}

motif_4 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (a1 - a2 - g * th + sqrt(4 * a1 * g * th + (a1 - a2 - g * th)^2)) /
    (2 * g)
  x2 <- (-a1 + a2 - g * th + sqrt(4 * a1 * g * th + (a1 - a2 - g * th)^2)) /
    (2 * g)
  
  return(c(x1, x2))
}

motif_5 <- function(a1, a2, g = 1, th = 5) {
  x1 <- (a1 * a2 - (g^2) * (th^2)) / (g * (a2 + g * th))
  x2 <- (a1 * a2 - (g^2) * (th^2)) / (g * (a1 + g * th))
  
  return(c(x1, x2))
}

motif_6 <- function(a1, a2, g = 1, th = 5) {
  x1 <- a1 / g
  x2 <- a2 / g
  
  return(c(x1, x2))
}

legacy_GtoP <- list(motif_1, motif_2, motif_3, motif_4, motif_5, motif_6)
J1 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- g*x1
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c(
    1/g,
    a2 * th / (g * (a1 + th)^2), 
    0, 
    a1/(g*(th + a1))
  ), nrow = 2))
}

J2 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- g*x1
  a2 <- (g*x2)/(1-(x1/(th+x1)))
  
  x1 <- a1/g
  x2 <- (a2*(1-((a1/g)/(th+((a1)/(g))))))/(g)
  
  
  return(matrix(c(
    1/g, 
    - a2 * th/(a1 + g * th)^2, 
    0, 
    (1-((a1/g)/(th+((a1)/(g)))))/(g)
  ), nrow = 2))
}

J3 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1)/(1-((x2)/(th+x2)))
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c((th*(a1 + 2*a2+g*th + sqrt(4 * a1 * a2 +(a1+g*th)^2)))/(2*(a2 + g * th)*sqrt(4 * a1 * a2 + (a1 + g * th)^2)), 
                  ((a1 + 2 * a2 + g * th)/sqrt(4 * a1 * a2 + (a1 + g * th)^2) -1)/(2*g), 
                  -(4 * a1^2 * th)/(sqrt(4 * a1 * a2 + (a1 + g * th)^2)*(-a1+g * th + sqrt(4 * a1 * a2 + (a1 + g * th)^2))^2), 
                  a1/(g * sqrt(4 * a1 * a2 + (a1 + g * th)^2))), nrow = 2))
}

J4 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(th)
  a2 <- (g*x2*(th+x1))/(th)
  
  return(matrix(c(
    (1+(a1 - a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (-1+(a1 - a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (-1+(-a1 + a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g), 
    (1+(-a1 + a2 + g * th)/(sqrt(4 * a1 * g * th + (-a1 + a2 + g * th)^2)))/(2*g)
  ), nrow = 2))
}

J5 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- (g*x1*(th+x2))/(x2)
  a2 <- (g*x2*(th+x1))/(x1)
  
  return(matrix(c(
    a2/(g*(a2 + th)), 
    (th*(a2 + th))/(g*(a1 + th)^2), 
    (th*(a1 + th))/(g*(a2 + th)^2), 
    a1/(g * (a1 + th))
  ), nrow = 2))
}

J6 <- function(x1, x2, g = 1, th = 5) {
  
  a1 <- x1*g
  a2 <- x2*g
  
  return(matrix(c(1/g, 0, 0, 1/g), nrow = 2))
}


legacy_S <- list(J1, J2, J3, J4, J5, J6)
motifPG_1 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2 * (th + x1)) / (x1)
  return(c(a1, a2))
}

motifPG_2 <- function(x1, x2, g = 1, th = 5) {
  a1 <- g * x1
  a2 <- (g * x2) / (1 - (x1 / (th + x1)))
  return(c(a1, a2))
}

motifPG_3 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1) / (1 - ((x2) / (th + x2)))
  a2 <- (g * x2 * (th + x1)) / (x1)
  
  return(c(a1, a2))
}

motifPG_4 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (th)
  a2 <- (g * x2 * (th + x1)) / (th)
  return(c(a1, a2))
}

motifPG_5 <- function(x1, x2, g = 1, th = 5) {
  a1 <- (g * x1 * (th + x2)) / (x2)
  a2 <- (g * x2 * (th + x1)) / (x1)
  return(c(a1, a2))
}

motifPG_6 <- function(x1, x2, g = 1, th = 5) {
  a1 <- x1 * g
  a2 <- x2 * g
  return(c(a1, a2))
}

motifPG <- function(k, x) list(motifPG_1, motifPG_2, motifPG_3, motifPG_4, motifPG_5, motifPG_6)[[k]](x[1], x[2])
