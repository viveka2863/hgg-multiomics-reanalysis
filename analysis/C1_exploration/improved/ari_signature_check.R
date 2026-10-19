ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab))); e <- sa * sb / c2(n)
  (s - e) / ((sa + sb) / 2 - e)
}
truth <- c(rep("Cortex", 3), rep("PD", 4), rep("TRK", 3))
cat("one TRK sample with the PD cluster:", ari(c(1,1,1,2,2,2,2,2,3,3), truth), "\n")
cat("one PD sample with the TRK cluster:", ari(c(1,1,1,2,2,2,3,3,3,3), truth), "\n")
cat("tumours only, one TRK with PD:", ari(c(1,1,1,1,1,2,2), c(rep("PD",4), rep("TRK",3))), "\n")
