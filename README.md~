Ts/Tv in HIV/SIV and EIV
========================


HIV/EIV
=======

Load data

```r
hiv_dat <- read.table("hiv/models_hiv.txt", head = T, as.is = T)


hiv_dat$file <- toupper(hiv_dat$file)

cases.pol <- grepl("POL|GAG", hiv_dat$file)
cases.env <- grepl("ENV", hiv_dat$file)
scales <- sapply(hiv_dat$file, function(x) strsplit(x, "[.]")[[1]][1])

tr_ts <- (hiv_dat$Q2 + hiv_dat$Q5)/(hiv_dat$Q1 + hiv_dat$Q3 + hiv_dat$Q4 + hiv_dat$Q6)
hiv_dat <- cbind(hiv_dat, tr_ts)
hiv_dat
```

```
##                   file   model   logLik         I k  shape    Q1     Q2
## 1        SC1.ENV.FASTA GTR+G+I -31714.8 0.0008915 4 0.8486 2.169  4.348
## 2        SC1.POL.FASTA GTR+G+I -10449.1 0.3299111 4 0.7365 1.835  9.776
## 3        SC2.ENV.FASTA GTR+G+I -53701.1 0.2613247 4 0.9429 1.647  3.956
## 4        SC2.POL.FASTA GTR+G+I -36038.9 0.3974838 4 0.9293 1.811  8.372
## 5      SC3.A.ENV.FASTA GTR+G+I -16197.8 0.2678737 4 0.6831 1.819  3.969
## 6      SC3.A.POL.FASTA GTR+G+I -21652.9 0.4429497 4 0.9012 1.321  4.630
## 7      SC3.B.ENV.FASTA GTR+G+I -24602.3 0.2481069 4 0.5546 2.376  4.934
## 8      SC3.B.POL.FASTA GTR+G+I -18446.1 0.5008064 4 0.8946 2.486 10.434
## 9        SC4.ENV.FASTA GTR+G+I -19212.7 0.3042801 4 0.6863 1.762  3.857
## 10       SC4.POL.FASTA GTR+G+I -13827.4 0.4418524 4 0.7916 2.561 10.367
## 11   SC5.ENV.IND.FASTA GTR+G+I  -9580.4 0.4331842 4 0.6409 2.070  4.340
## 12    SC5.ENV.SA.FASTA GTR+G+I -27580.2 0.4952168 4 0.8050 1.773  4.208
## 13    SC5.POL.SA.FASTA GTR+G+I  -7536.5 0.5081877 4 0.9163 3.270 10.672
## 14 SC6.ENV.PAT85.FASTA GTR+G+I  -5809.7 0.7990771 4 0.7168 2.222  4.417
## 15 SC6.GAG.PAT21.FASTA GTR+G+I   -954.2 0.6594123 4 1.0170 1.632 12.469
## 16 SC6.POL.PAT53.FASTA GTR+G+I  -1663.2 0.2380274 4 1.0030 2.329 12.851
## 17 SC6.POL.PAT66.FASTA GTR+G+I  -1801.8 0.7794093 4 0.8148 2.494 21.735
##           Q3        Q4     Q5 Q6    bf1    bf2    bf3    bf4 meanKaks
## 1  1.017e+00 1.087e+00  5.126  1 0.3681 0.1835 0.2159 0.2325   1.0770
## 2  1.419e+00 1.338e+00 15.962  1 0.4229 0.1793 0.2004 0.1974   1.2112
## 3  6.301e-01 7.662e-01  4.430  1 0.3813 0.1767 0.2138 0.2282   3.2003
## 4  8.148e-01 8.065e-01 11.698  1 0.4125 0.1615 0.2072 0.2188   5.5589
## 5  6.213e-01 7.121e-01  4.156  1 0.3560 0.1764 0.2267 0.2409   0.5956
## 6  5.339e-01 4.876e-01  6.620  1 0.3795 0.1723 0.2335 0.2147   0.2521
## 7  9.414e-01 1.163e+00  4.941  1 0.3660 0.1643 0.2312 0.2384   1.6767
## 8  9.534e-01 7.703e-01 12.386  1 0.3911 0.1688 0.2251 0.2150   1.1776
## 9  4.235e-01 6.496e-01  4.091  1 0.3591 0.1780 0.2211 0.2418   1.1913
## 10 9.717e-01 7.672e-01 11.781  1 0.3878 0.1666 0.2273 0.2184   3.4420
## 11 6.194e-01 7.871e-01  3.966  1 0.3530 0.1705 0.2350 0.2416   0.4716
## 12 7.731e-01 5.580e-01  5.075  1 0.3618 0.1771 0.2393 0.2217   0.4110
## 13 9.467e-01 7.899e-01 15.789  1 0.3910 0.1681 0.2273 0.2136   0.3007
## 14 5.872e-01 9.174e-01  4.064  1 0.3508 0.1730 0.2359 0.2403   0.2169
## 15 3.262e+00 2.106e-05 17.429  1 0.3884 0.1923 0.2451 0.1742   0.2454
## 16 5.373e-01 9.859e-06  5.173  1 0.3940 0.1626 0.2072 0.2363   0.0000
## 17 6.218e-06 3.398e-05 19.803  1 0.3875 0.1653 0.2095 0.2376   0.0750
##      varKaks minTime maxTime  tr_ts
## 1  2.869e-01    1959    2009  1.797
## 2  9.026e-01    1960    2007  4.603
## 3  4.573e+05    1983    2009  2.074
## 4  2.294e+03    1983    2008  4.528
## 5  1.575e+00    1991    2008  1.957
## 6  3.466e-01    1992    2005  3.366
## 7  2.769e+02    1983    2004  1.802
## 8  1.387e+02    1983    2009  4.381
## 9  6.102e+01    1991    2008  2.072
## 10 4.937e+02    1985    2006  4.179
## 11 8.520e-01    1999    2001  1.856
## 12 5.030e-01    1992    2007  2.262
## 13 3.191e-01    1992    2007  4.405
## 14 2.400e-01    2002    2005  1.794
## 15 1.542e-01    2004    2008  5.072
## 16 0.000e+00    1999    2000  4.662
## 17 5.987e-02    1999    2002 11.887
```



Run plots with log10 transformation


```r

tr_ts_scale <- cbind(factor(scales), tr_ts)

par(mar = c(5, 5, 4, 4))
### Tr Ts ratio plot
plot(jitter(as.numeric(gsub("SC", "", scales))), log10(tr_ts), col = c(rgb(1, 
    0, 0, 0.8), rgb(0, 0, 1, 0.8))[cases.pol + 1], xlim = c(6.8, 0), pch = 20, 
    axes = F, ann = T, type = "p", ylab = expression(log[10](italic(Ts/Tv))), 
    xlab = expression(paste(bold("Fig1. "), (italic(Ts/Tv)), " by sampling level in HIV and SIV data")), 
    ylim = c(0.1, 1.5), cex = 3)

axis(1, lab = F)
axis(2)
text(axTicks(1), -0.02, labels = c("", "Transmission\nchain", "Between\nhosts", 
    "Host\npopulation", "Between\npopulations", "Global\ndiversity", "Interspecies", 
    ""), xpd = T, cex = 0.7)
legend(1, 2, legend = c(expression(italic(ENV)), expression(italic(POL))), text.col = c(rgb(1, 
    0, 0, 0.8), rgb(0, 0, 1, 0.8)), bty = "n")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 





```r
### KA Ks ratio plot
plot(jitter(as.numeric(gsub("SC", "", scales))), log10(hiv_dat$meanKaks), col = c(rgb(1, 
    0, 0, 0.8), rgb(0, 0, 1, 0.8))[cases.pol + 1], xlim = c(6.8, 0), pch = 20, 
    axes = F, ann = T, type = "p", ylab = expression(log[10](italic(dN/dS))), 
    xlab = expression(paste(bold("Fig2. "), italic(dN/dS), " by sampling level in HIV and SIV data")), 
    cex = 3)

axis(1, lab = F)
axis(2)
text(axTicks(1), -1.3, labels = c("", "Transmission\nchain", "Between\nhosts", 
    "Host\npopulation", "Between\npopulations", "Global\ndiversity", "Interspecies", 
    ""), xpd = T, cex = 0.7)
legend(1, -0.5, legend = c(expression(italic(ENV)), expression(italic(POL))), 
    text.col = c(rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)), bty = "n")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 






```r
## TsTv dNdS plot
par(mar = c(6, 5, 5, 5))
plot(log10(hiv_dat$tr_ts), log10(hiv_dat$meanKaks), col = c(rgb(1, 0, 0, 0.8), 
    rgb(0, 0, 1, 0.8))[cases.pol + 1], , pch = 20, axes = F, ann = T, type = "p", 
    ylab = expression(log[10](italic(dN/dS))), xlab = expression(log[10](italic(Ts/Tv))), 
    ylim = c(-1.1, 0.6), xlim = c(0.1, 1.2), cex = 3)

axis(1)
axis(2)
mtext(expression(paste(bold("Fig 3. "), italic(dN/dS), " vs.", italic(Ts/Tv))), 
    side = 1, line = 4.5)

legend(1, -0.5, legend = c(expression(italic(ENV)), expression(italic(POL))), 
    text.col = c(rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)), bty = "n")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


EIV
====

Load EIV data


```r
eiv_dat <- read.table("EIV/models_hiv.txt", head = T)


tr_ts_eiv <- (eiv_dat$Q2 + eiv_dat$Q5)/(eiv_dat$Q1 + eiv_dat$Q3 + eiv_dat$Q4 + 
    eiv_dat$Q6)
eiv_dat <- cbind(eiv_dat, tr_ts_eiv)
eiv_dat <- eiv_dat[c(3, 2, 1), ]

print(eiv_dat)
```

```
##                       file   model logLik         I k   shape Q1     Q2 Q3
## 3    EIV.HA.Outbreak.fasta GTR+G+I  -1265 5.410e-05 4 499.994  1  8.128  1
## 2 EIV.HA.Intrahost-2.fasta GTR+G+I  -1578 5.684e-05 4 499.994  1  8.613  1
## 1    EIV.HA.Global-3.fasta GTR+G+I  -5710 3.161e-01 4   1.218  1 14.359  1
##   Q4     Q5 Q6    bf1    bf2    bf3    bf4 meanKaks varKaks minTime
## 3  1  8.128  1 0.3582 0.2015 0.2066 0.2336   0.0500 0.05000      NA
## 2  1  8.613  1 0.3536 0.2033 0.2115 0.2316   0.0716 0.04172    2001
## 1  1 14.359  1 0.3640 0.1879 0.2077 0.2403   0.1214 0.10492    1963
##   maxTime tr_ts_eiv
## 3      NA     4.064
## 2    2012     4.306
## 1    2008     7.179
```



Run similar plots to those for the previous data set



```r
par(mar = c(5, 5, 4, 4))
### Tr Ts ratio plot

plot(c(1, 2, 3), log10(eiv_dat$tr_ts), col = c(rgb(1, 0, 0, 0.8)), type = "p", 
    pch = 20, axes = F, ann = T, ylab = expression(log[10](italic(Ts/Tv))), 
    xlab = expression(paste(bold("Fig4. "), (italic(Ts/Tv)), " by sampling level in the EIV data")), 
    cex = 3)
axis(1, lab = F)
text(axTicks(1), 0.58, labels = c("Outbreak", "", "Intrahost", "", "Global"), 
    xpd = T, cex = 0.7)
axis(2)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 




```r
### Ka Ks ratio plot

plot(c(1, 2, 3), log10(eiv_dat$meanKaks), col = c(rgb(1, 0, 0, 0.8)), pch = 20, 
    axes = F, ann = T, type = "p", ylab = expression(log[10](italic(dN/dS))), 
    xlab = expression(paste(bold("Fig5. "), italic(dN/dS), " by sampling level in the EIV data")), 
    cex = 3)
axis(1, lab = F)
axis(2)
text(axTicks(1), -1.35, labels = c("Outbreak", "", "Intrahost", "", "Global"), 
    xpd = T, cex = 0.7)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 




```r
### TsTv dNdS plot

par(mar = c(6, 5, 5, 5))
plot(log10(eiv_dat$tr_ts), log10(eiv_dat$meanKaks), col = rgb(1, 0, 0, 0.8), 
    pch = 20, axes = F, ann = T, type = "p", ylab = expression(log[10](italic(dN/dS))), 
    xlab = expression(log[10](italic(Ts/Tv))), cex = 3)

axis(1)
axis(2)
mtext(expression(paste(bold("Fig 6. "), italic(dN/dS), " vs.", italic(Ts/Tv))), 
    side = 1, line = 4.5)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 



