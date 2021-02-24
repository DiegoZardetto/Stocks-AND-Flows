##############################################################
# One prototype of complete Demographic Balancing procedure. #
##############################################################

# S T A R T #

# Timestamp
Tstart <- print(format(Sys.time(), "%A %d %B %Y | %X"))

# Source the Balancing code
source("G:\\BalanceR\\Engine\\Balancer.R")

# Source the DBE layer code
source("G:\\BalanceR\\DBE\\StandardInputCode.R")

#############################################################################
# PHASE 1: Read the input data (in Standard Format) from external files and #
#          prepare them for the Balancing process.                          #
#############################################################################

# Directory with standard input data
setwd("G:\\BalanceR\\DBE\\Standard_Input")

# Initial population counts
P1 <- read.P(file = "P1.csv")

# Final population counts
P2 <- read.BDN(file = "P2.csv", STOCK = P1)

# Births
B <- read.BDN(file = "B.csv", STOCK = P1)

# Deaths
D <- read.BDN(file = "D.csv", STOCK = P1)

# Deduce Natural Increase N = B - D
N <- B
N$N <- B$N - D$N

# (Generalized) Migration Flows
F <- read.F(file = "F.csv", STOCK = P1)

# Structural constraints on (Generalized) Migration Flows
Cmask <- mkCmask(F)

# Test structural constraints on the (Generalized) Migration Flows Matrix
mkCheck(F, Cmask)

# Deduce the Net (Generalized) Migration Flows Matrix
M <- t(F) - F

# Settle the 'Abroad' values
# A) The TRADITIONAL approach: settle P1 to *0* (as it already is after the input phase)
#    and derive P2 accordingly from the DBE 
P1 <- settle.AbroadP1(P1, setval = 0)
P2 <- settle.AbroadP2(P2, P1, N, M)

# Free some memory (hopefully)
gc()

# Timestamp
Tread <- print(format(Sys.time(), "%A %d %B %Y | %X"))


#############################################################################
# PHASE 2: Transform input data according to the input format required by   #
#          BalanceR, and solve the demographic Balancing problem.           #
#############################################################################

# Build the objects to be balanced
o <- makeObj(P1, P2, B, D, N, F, M)

# Build the variance structure of the objects to be balanced.
# Possible models:
# 1) POISSON
# 2) POISSON-SKELLAM
# 3) HOMOSKEDASTIC

# NOTE: If needed (e.g. if varmodel = "poisson-skellam"), choose
#       how to treat 0 raw counts (by default, M.N.adj.zeros = FALSE)
# ov <- makeObjV(o, Cmask, varmodel = "poisson")
# ov <- makeObjV(o, Cmask, varmodel = "poisson-skellam", M.N.adj.zeros = FALSE)
ov <- makeObjV(o, Cmask, varmodel = "poisson-skellam", M.N.adj.zeros = TRUE)
# ov <- makeObjV(o, Cmask, varmodel = "poisson-skellam", M.N.adj.zeros = TRUE, vD = 1, vN = 1)
# ov <- makeObjV(o, Cmask, varmodel = "poisson-skellam", M.N.adj.zeros = FALSE, vP2 = 0, vB = 1, vD = 1, vN = 1)
# ov <- makeObjV(o, Cmask, varmodel = "poisson-skellam", M.N.adj.zeros = FALSE, vP2 = 0)

# Read the balancing constraints
bc <- read.BalConstr(file = "G:\\BalanceR\\DBE\\DBEbalconstr.txt")

# Set the objects (w.r.t. constraints)
bo <- set.BalObj(o, bc)

# Set the variances (w.r.t. constraints)
bov <- set.BalObj(ov, bc)

# Compute balancing results
boq <- Balance(bc, bo, bov, verbose = TRUE)

# Free some memory (hopefully)
gc()

# Timestamp
Tbal <- print(format(Sys.time(), "%A %d %B %Y | %X"))


#############################################################################
# PHASE 3: Write the balanced output data (in Standard Format) to external  #
#          files.                                                           #
#############################################################################

# Directory to write standard output data
# outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\POISSON"
# outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\POISSON_SKELLAM"
outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\POISSON_SKELLAM_AdjZeros"
# outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\DEATHS_POISSON_SKELLAM_AdjZeros"
# outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\BIRTHS_DEATHS_NO_P2_POISSON_SKELLAM"
# outdir <- "G:\\BalanceR\\DBE\\Standard_Output\\NO_P2_POISSON_SKELLAM"

# Move the working directory there (in case IO operations are needed)
setwd(outdir)

# Write ALL balanced objects in standard format
Q <- write.ALL(boq, STOCK = P1, FLOWS = F, path = outdir)
# Q

# Compare balanced and original objects and write to external files
Delta <- compare.ALL(boq, bo, STOCK = P1, FLOWS = F, path = outdir)
# Delta

# Free some memory (hopefully)
gc()

# Timestamp
Twrite <- print(format(Sys.time(), "%A %d %B %Y | %X"))


###########################################################################
# PHASE 4: Run diagnostics on the obtained output.                        #
#          Typically: (1) check residual discrepancies in atomic DBEs.    #
#                     (2) check for negative balanced estimates of counts #
#                         that should be positive.                        #
#                     (3) inspect the relevant atomic DBEs to trace the   #
#                         likely origin of "pathologies" (1) and (2).     #
###########################################################################

# Assess discrepancies in atomic DBEs for raw estimates
Errors.raw <- getDBEerrors(bo, P1, P2, N, M, tol = 1E-7)
# Have a look
# Errors.raw

# Assess possible residual discrepancies in atomic DBEs for balanced estimates
Errors.bal <- getDBEerrors(boq, P1, P2, N, M, tol = 1E-7)
# Have a look
# Errors.bal

# Use function MicroMonitor() to inspect the DBEs that are still not satisfied
# mm <- MicroMonitor(Delta, i = c(Sex1 = NA, Nation1 = NA, AgeCl1 = NA, Region1 = NA))

# Check for negative balanced estimates of counts that should be positive
Negatives <- checkNegativeEstimates(Delta)
# Have a look
Negatives

# Use function MicroMonitor() to inspect the DBEs related to negative estimates
# mm <- MicroMonitor(Delta, i = c(Sex1 = NA, Nation1 = NA, AgeCl1 = NA, Region1 = NA))

# Free some memory (hopefully)
gc()

# Timestamp
Tdiagn <- print(format(Sys.time(), "%A %d %B %Y | %X"))

# Save the session history in the standard output directory 
savehistory(file = "CompleteLog.Rhistory")

# Save the current workspace in the standard output directory 
save.image(file = "Output.RData")

# Timestamp
Tend <- print(format(Sys.time(), "%A %d %B %Y | %X"))

# E N D #
