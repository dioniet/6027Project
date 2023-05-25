#samples of fish collected in the Bay of Seine (Galgani et al. 1991).
#https://www.sciencedirect.com/science/article/pii/0025326X9190403F

Y0 <- matrix(
  c(9100, 98, 38, 81, 1.7, 15, 17,
    2000, 26, 11, 12, 1.2, 10, 9.2,
    6700, 68, 27, 59, 1.3, 10, 7.1,
    3600, 64, 24, 11, 2.9, 22, 2.3,
    2300, 44, 17, 13, 2.0, 16, 11,
    3200, 51, 12, 20, 2.4, 22, 50,
    2600, 41, 15, 7.6,1.1, 8.0,5.3,
    13000,190,37, 23, 0.9, 8.7,14,
    5800, 71, 29, 56, 1.5, 13, 14,
    2800, 45, 14, 26, 1.9, 10, 12),
  nrow = 10, byrow = T, dimnames = list(rownames =
                                          c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                                        colnames =
                                          c("C1:PCB", "C2: DDE", "C3: DDD", "C4: DDT", "C5: α-HCH", "C6:γ-HCH","C7: PAH")))
Y0

