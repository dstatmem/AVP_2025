library(readr)
CA_prec <- read_delim("C:/Users/admin/Dropbox/MSc Maria/Maria/Field/Data_files/CA_prec.csv", 
                        +     ";", escape_double = FALSE, trim_ws = TRUE)

plot(lm(CA_prec$Ca~CA_prec$`Alk ueq`))
summary(lm(CA_prec$Ca~CA_prec$`Alk ueq`))
summary(lm(CA_prec$`Alk ueq`~CA_prec$Ca))

plot(CA_prec$`Alk ueq` ~ CA_prec$Ca, ylab = "Alkalinity [?? equivalent/L]", xlab = "Caclium [umol/L]",  pch = 16, col = "darkblue")
lines(127.22657 + 2.17576 * CA_prec$Ca ~ CA_prec$Ca , col = "red")
lines(127.22657 + (2.17576 + 0.1619156)* CA_prec$Ca ~ CA_prec$Ca , col = "red", lty = 3)
lines(127.22657 + (2.17576 - 0.1619156)* CA_prec$Ca ~ CA_prec$Ca , col = "red", lty = 3)
legend("topleft", 
       inset = c(0.05,0.05), 
       cex = 1.1, 
       bty = "n", 
       legend = c("Alkalinity ??eq", "Regresion line", "95% confiedence interval"), 
       text.col = c("black", "black", "black"),
       col = c("darkblue", "red", "red"), 
       lwd=1, lty=c(NA,1,3), pch=c(16,NA,NA))
