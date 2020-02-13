library(readr)

a <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Deposit 1.tsv")

a1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Deposit 1.tsv")

test1<-a == a1


b <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Deposit 7b.tsv")

b1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Deposit 7b.tsv")

test2 <-b == b1


c <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Deposit 13.tsv")

c1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Deposit 13.tsv")

test3 <-c == c1

d <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Deposit 14.tsv")

d1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Deposit 14.tsv")

test <- d == d1



e <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Misc Deposit 1.tsv")

e1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Misc Deposit 1.tsv")

test4<-e == e1


f <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Misc Deposit 7b.tsv")

f1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Misc Deposit 7b.tsv")

test5<-f == f1


g <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Misc Deposit 13.tsv")

g1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Misc Deposit 13.tsv")

test6 <-g == g1


h <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Misc Deposit 14.tsv")

h1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Misc Deposit 14.tsv")

test7 <-h == h1

i <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/Old_versions/RLB_P23_Small_Mammals - Misc Hancock.tsv")

i1 <- read_tsv("data/original_google_data/GoogleDriveExports-mammals/RLB_P23_Small_Mammals - Misc Hancock.tsv")

test8 <-i == i1



