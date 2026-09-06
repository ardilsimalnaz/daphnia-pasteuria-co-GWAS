K <- 2:10

daphnia_delta <- c(
  0.00683,
  0.00204,
  0.00330,
  0.00046,
  0.00155,
 -0.00113,
 -0.00194,
  0.00014,
  0.00008
)

pasteuria_delta <- c(
  0.44549,
  0.21354,
  0.05356,
  0.02697,
  0.00940,
  0.00412,
 -0.00949,
  0.01532,
 -0.00706
)

png("ADMIXTURE_deltaCV.png",
    width=1600, height=700, res=150)

par(mfrow=c(1,2),
    mar=c(5,5,3,1))

plot(K, daphnia_delta,
     type="b",
     pch=19,
     xlab="K",
     ylab="Change in CV error",
     main="Daphnia magna",
     xaxt="n")

axis(1, at=K)
abline(h=0, lty=2)

plot(K, pasteuria_delta,
     type="b",
     pch=19,
     xlab="K",
     ylab="Change in CV error",
     main="Pasteuria ramosa",
     xaxt="n")

axis(1, at=K)
abline(h=0, lty=2)

dev.off()
