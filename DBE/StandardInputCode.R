################################
# S T A N D A R D    I N P U T #
################################

`read.P` <- function(file = if (interactive()) choose.files(caption = "Select POPULATION STOCKS file", multi = FALSE, filters = Filters["All", ]),
                     sep = ";", dec = '.', ...) {
##########################################################################
# Read input POPULATION COUNTS from external files and prepare them.     #
#                                                                        #
# NOTE: Input data are assumed to comply with the STANDARD INPUT FORMAT. #
#                                                                        #
# NOTE: The prepared output has a STANDARD ORDERING of rows and columns. #
##########################################################################
  # Set the MAXIMAL schema for POPULATION STOCKS
  stockClasses <- c("Sex" = "factor", "Nation" = "factor", "AgeCl" = "factor", "Region" = "factor", "N" = "numeric")
  stockSchema <- names(stockClasses)[stockClasses == "factor"]

  # First read the data, then prepare it
  P <- suppressWarnings(read.table(file, sep = sep, dec = dec, head = TRUE, colClasses = stockClasses, ...))

  # If there are unknown variables, raise an error
  unknown.vars <- names(P)[!(names(P) %in% names(stockClasses))]
  if (length(unknown.vars) > 0) {
     stop("Unknown variables in POPULATION STOCKS file: ", paste(unknown.vars, collapse = ", "), ".\n")
    }

  # Reduce the stockSchema to the REALIZED one
  stockSchema <- stockSchema[stockSchema %in% names(P)]

  # IF NEEDED, add 'abroad' level 'E' to Region (it will be the LAST level, as
  # according to the STANDARD INPUT FORMAT Region's levels are integers 1, ..., k)
  if ("Region" %in% stockSchema) {
     levels(P$Region) <- c(levels(P$Region), "E")
    }

  # IF NEEDED, add 'newborn' level 'N' to AgeCl, and set it as the FIRST level
  if ("AgeCl" %in% stockSchema) {
     levels(P$AgeCl) <- c(levels(P$AgeCl), "N")
     P$AgeCl <- relevel(P$AgeCl, "N")
    }

  # Tidy up the levels of stockSchema's factors
  for (fact in stockSchema) {
     P[[fact]] <- factor(P[[fact]], levels = tidyLevels(P[[fact]]))
    }

  # Build ALL the rows of the stockSchema, including empty levels
  stockSchemaRows <- expand.grid(lapply(P[, stockSchema, drop = FALSE], levels))

  # Merge P with stockSchemaRows (note that this does NOT preserve FACTOR ordering, in general)
  # REMARK: As (i)  stockSchemaRows and P have only common factors, and
  #            (ii) stockSchemaRows comes first in merge(), THE COLUMNS OF P WILL APPEAR
  #         IN THE SAME ORDER AS THEY APPEAR IN stockSchema
  P <- merge(stockSchemaRows, P, all = TRUE)

  # Order the rows of P so that the RIGHTMOST variable varies FASTER
  P <- P[do.call(order, P[, stockSchema, drop = FALSE]), ]

  # IF NEEDED, put to 0 all the counts when AgeCl1 = 'N' or Region1 = 'E'
  P$N[is.na(P$N)] <- 0

  # Return the prepared object
  return(P)
}


`read.BDN` <- function(file = if (interactive()) choose.files(caption = "Select BIRTHS / DEATHS / NATURAL INCREASE file", multi = FALSE, filters = Filters["All", ]),
                       sep = ";", dec = '.', STOCK, ...) {
#############################################################################
# Read input BIRTHS / DEATHS / NATURAL INCREASE figures from external files #
# and prepare them.                                                         #
#                                                                           #
# NOTE: Input data are assumed to comply with the STANDARD INPUT FORMAT.    #
#                                                                           #
# NOTE: The prepared output has a STANDARD ORDERING of rows and columns.    #
#                                                                           #
# NOTE: STOCK must be an already prepared POPULATION COUNTS object.         #
#       It serves the purpose to guarantee stocks-flows FORMAT consistency. #
#                                                                           #
# REMARK: This function may *and should* be used to read T2 (T1) POPULATION #
#         COUNTS while using for argument STOCK already prepared T1 (T2)    #
#         POPULATION COUNTS. This would ensure P1 <-> P2 congruence!        #
#############################################################################
  # Set the MAXIMAL schema for POPULATION STOCKS
  stockClasses <- c("Sex" = "factor", "Nation" = "factor", "AgeCl" = "factor", "Region" = "factor", "N" = "numeric")
  stockSchema <- names(stockClasses)[stockClasses == "factor"]

  # First read the data, then prepare it
  BDN <- suppressWarnings(read.table(file, sep = sep, dec = dec, head = TRUE, colClasses = stockClasses, ...))

  # If there are unknown variables, raise an error
  unknown.vars <- names(BDN)[!(names(BDN) %in% names(stockClasses))]
  if (length(unknown.vars) > 0) {
     stop("Unknown variables in BIRTHS / DEATHS / NATURAL INCREASE file: ", paste(unknown.vars, collapse = ", "), ".\n")
    }

  # Reduce the stockSchema to the REALIZED one
  stockSchema <- stockSchema[stockSchema %in% names(BDN)]

  # IF NEEDED, add 'abroad' level 'E' to Region (it will be the LAST level, as
  # according to the STANDARD INPUT FORMAT Region's levels are integers 1, ..., k)
  if ("Region" %in% stockSchema) {
     levels(BDN$Region) <- c(levels(BDN$Region), "E")
    }

  # IF NEEDED, add 'newborn' level 'N' to AgeCl, and set it as the FIRST level
  if ("AgeCl" %in% stockSchema) {
     levels(BDN$AgeCl) <- c(levels(BDN$AgeCl), "N")
     BDN$AgeCl <- relevel(BDN$AgeCl, "N")
    }

#  # Tidy up the levels of stockSchema's factors
#  for (fact in stockSchema) {
#     BDN[[fact]] <- factor(BDN[[fact]], levels = tidyLevels(BDN[[fact]]))
#    }

  # Tidy up the levels of stockSchema's factors. Try ensure the congruence
  # of BIRTHS / DEATHS / NATURAL INCREASE matrices with STOCKS matrices
  # NOTE: Code below implies that:
  #       (i)  factor combinations that are possibly *absent* in input will
  #            appear as zero frequency cells *present* in output object BDN,
  #            which is *correct*.
  #       (ii) factor combinations that are *present* in input but possibly *absent*
  #            in STOCK would appear as *extra rows* in output object BDN,
  #            which triggers a *preventive informative error message*.
  for (fact in stockSchema) {
     BDN[[fact]] <- factor(BDN[[fact]], levels = levels(STOCK[[fact]]))
     if (anyNA(BDN[[fact]])) {
         stop("Found NA values in input BIRTHS / DEATHS / NATURAL INCREASE file's factor: ", fact,
              "\n  NOTE: This might indicate that the input file involves levels which are absent in STOCK!")
        }
    }

  # Get from STOCK the stockSchemaRows, i.e. ALL the needed rows including empty levels
  stockSchemaRows <- STOCK[, stockSchema, drop = FALSE]

  # Merge BDN with stockSchemaRows (note that this does NOT preserve FACTOR ordering, in general)
  # REMARK: As (i)  stockSchemaRows and BDN have only common factors, and
  #            (ii) stockSchemaRows comes first in merge(), THE COLUMNS OF BDN WILL APPEAR
  #         IN THE SAME ORDER AS THEY APPEAR IN stockSchema
  BDN <- merge(stockSchemaRows, BDN, all = TRUE)

  # Order the rows of BDN so that the RIGHTMOST variable varies FASTER
  BDN <- BDN[do.call(order, BDN[, stockSchema, drop = FALSE]), ]

  # IF NEEDED, put to 0 all the counts when AgeCl1 = 'N' or Region1 = 'E'
  BDN$N[is.na(BDN$N)] <- 0

  # Check for possible dimension mismatches between BDN output and input STOCK matrix
  # NOTE: This error should *never* pop-up, as it should have already been catched above,
  #       when tidying up flowSchema's factors
  if ( NROW(STOCK) != NROW(BDN) ) {
     stop("Number of rows of BIRTHS / DEATHS / NATURAL INCREASE matrix (", NROW(BDN), ") does not agree with STOCK's one (", NROW(STOCK), ")!")
    }

  # Return the prepared object
  return(BDN)
}


`read.F` <- function(file = if (interactive()) choose.files(caption = "Select MIGRATION FLOWS file", multi = FALSE, filters = Filters["All", ]),
                     sep = ";", dec = '.', STOCK, ...) {
#############################################################################
# Read input MIGRATION FLOWS from external files and prepare them.          #
#                                                                           #
# NOTE: Input data are assumed to comply with the STANDARD INPUT FORMAT.    #
#                                                                           #
# NOTE: The prepared output has a STANDARD ORDERING of rows and columns.    #
#                                                                           #
# NOTE: STOCK must be an already prepared POPULATION COUNTS object.         #
#       It serves the purpose to guarantee stocks-flows FORMAT consistency. #
#############################################################################
  # Set the MAXIMAL schema for MIGRATION FLOWS
  flowClasses <- c("Sex1" = "factor", "Nation1" = "factor", "AgeCl1" = "factor", "Region1" = "factor",
                   "Sex2" = "factor", "Nation2" = "factor", "AgeCl2" = "factor", "Region2" = "factor", "N" = "numeric")  
  flowSchema <- names(flowClasses)[flowClasses == "factor"]

  # First read the data, then prepare it
  Fin <- suppressWarnings(read.table(file, sep = sep, dec = dec, head = TRUE, colClasses = flowClasses, ...))

  # If there are unknown variables, raise an error
  unknown.vars <- names(Fin)[!(names(Fin) %in% names(flowClasses))]
  if (length(unknown.vars) > 0) {
     stop("Unknown variables in MIGRATION FLOWS file: ", paste(unknown.vars, collapse = ", "), ".\n")
    }

  # Reduce the flowSchema to the REALIZED one
  flowSchema <- flowSchema[flowSchema %in% names(Fin)]

  # IF NEEDED, add 'abroad' level 'E' to Region (it will be the LAST level, as
  # according to the STANDARD INPUT FORMAT Region's levels are integers 1, ..., k)
  if ("Region1" %in% flowSchema) {
     levels(Fin$Region1) <- c(levels(Fin$Region1), "E")
     levels(Fin$Region2) <- c(levels(Fin$Region2), "E")
    }

  # IF NEEDED, add 'newborn' level 'N' to AgeCl, and set it as the FIRST level
  if ("AgeCl1" %in% flowSchema) {
     levels(Fin$AgeCl1) <- c(levels(Fin$AgeCl1), "N")
     Fin$AgeCl1 <- relevel(Fin$AgeCl1, "N")
     levels(Fin$AgeCl2) <- c(levels(Fin$AgeCl2), "N")
     Fin$AgeCl2 <- relevel(Fin$AgeCl2, "N")
    }

#  # Tidy up the levels of flowSchema's factors
#  for (fact in flowSchema) {
#     Fin[[fact]] <- factor(Fin[[fact]], levels = tidyLevels(Fin[[fact]]))
#    }

  # Tidy up the levels of flowSchema's factors. Try ensure the congruence
  # of FLOWS matrices with STOCKS matrices
  # NOTE: Code below implies that:
  #       (i)  factor combinations that are possibly *absent* in Fin will appear
  #            as zero frequency cells *present* in output ftable object F, which
  #            is *correct*.
  #       (ii) factor combinations that are *present* in Fin but possibly *absent*
  #            in STOCK would appear as *extra rows/cols* in output ftable object F,
  #            which triggers a *preventive informative error message*.
  for (fact in flowSchema) {
     Fin[[fact]] <- factor(Fin[[fact]], levels = levels(STOCK[[substr(fact, 1, nchar(fact) - 1)]]))
     if (anyNA(Fin[[fact]])) {
         stop("Found NA values in input FLOWS file's factor: ", fact,
              "\n  NOTE: This might indicate that FLOWS file involves levels which are absent in STOCK!")
        }
    }

  # Expand the migration flows table using the count column 'N' in order to later
  # use ftable
  FFin <- Fin[rep(1:nrow(Fin), Fin$N), names(Fin) != "N"]

  # Compute the flat contingency table
  # REMARK: ftable() creates combinations of factor levels in such a way that the levels
  #         of the left-most variable vary the slowest. THIS MATCHES THE STANDARD ORDER
  #         created by functions read.P() and read.BDN().
  F <- ftable(FFin, row.vars = names(FFin)[endsWith(names(FFin), "1")],
                    col.vars = names(FFin)[endsWith(names(FFin), "2")], exclude = NULL)

  # Check for possible dimension mismatches between FLOWS output ftable and input STOCK matrix
  # NOTE: This error should *never* pop-up, as it should have already been catched above,
  #       when tidying up flowSchema's factors
  if ( ( NROW(STOCK) != NROW(F) ) ||  ( NROW(STOCK) != NCOL(F) ) ) {
     stop("FLOW matrix dimension (", NROW(F)," X ", NCOL(F), ") does not agree with STOCK's rows (", NROW(STOCK), ")!")
    }

  # Return the prepared object
  return(F)
}


`tidyLevels` <- function(factor) {
###########################################################################
# Arrange the levels of any classification variable belonging to the      #
# STANDARD INPUT (be it in stocks or flows) in STANDARD ORDER.            #
#                                                                         #
# The STANDARD ORDER of the levels is as follows:                         #
# 1) *Integer-like* levels, i.e. levels that can be read as INTEGER must  #
#    be *sorted* in ascending order.                                      #
#                                                                         #
# 2) *Character-like* levels, i.e. levels that can *NOT* be read as       #
#    INTEGER, can only be either:                                         #
#    2.1) "N" (for the newborn class in AgeCl) -> must be the FIRST level #
#    2.2) "E" (for the abroad class of Region) -> must be the LAST level  #
#                                                                         #
# 3) When *integer-like* and *character-like* levels are simultaneously   #
#    present, the standard order is obtained by applying rules 1) and 2). #
###########################################################################
  if (!is.factor(factor)) stop("Input variable must be a factor!")

  old.levels <- levels(factor)
  coerced.levels <- suppressWarnings(as.integer(old.levels))
  int.levels <- sort(coerced.levels[!is.na(coerced.levels)])
  char.levels <- old.levels[is.na(coerced.levels)]

  if (length(char.levels) > 1) stop("Found more than one character-like level in factor!")

  if (length(char.levels) == 0) {
     new.levels <- int.levels
    }

  if (length(char.levels) == 1) {
     if (!(char.levels %in% c("E", "N"))) stop("Found illicit character-like level in factor! (\"", char.levels, "\")")

     if (char.levels == "E") {
         # Thus factor must be Region: abroad level "E" must be the last
         new.levels <- c(int.levels, "E")
        }
     if (char.levels == "N") {
         # Thus factor must be AgeCl: newborn level "N" must be the first
         new.levels <- c("N", int.levels)
        }
    }

  if (length(new.levels) != length(old.levels)) stop("Length of tidy and old levels differ!")

  return(new.levels)
}


`mkCmask` <- function(F, print = FALSE) {
################################################################################
# Given a Migration Flows Matrix 'F', build its STRUCTURAL CONSTRAINTS mask,   #
# namely a flat contingency table of the same dimensions whose cell values     #
# are:                                                                         #
# A) 0 IFF the cell is FEASIBLE                                                #
# B) 1 IFF the cell is UNFEASIBLE                                              #
#                                                                              #
# NOTE: The counts of F are NEVER USED, ONLY the STRUCTURE of F MATTERS!       #
################################################################################
  # Get F classification variables
  row.vars <- names(attr(F, "row.vars"))
  col.vars <- names(attr(F, "col.vars"))
  vars <- c(row.vars, col.vars)

  # Convert F from ftable format to data.frame
  Cmask <- as.data.frame(F)

  # Initialize 'Freq' values to 2, then put to 1 those of UNFEASIBLE cells
  Cmask$Freq <- 2

  ############################
  ## STRUCTURAL CONSTRAINTS ##
  ############################
  ## C1) Sex CANNOT CHANGE
         ## NOTE: Even if it happens, the cases are so few that they WOULD
         ##       VERY LIKELY GENERATE ISSUES (e.g. negative balanced counts)
   if ("Sex1" %in% vars) {
         Cmask[ Cmask$Sex1 != Cmask$Sex2, "Freq"] <- 1
    }

  ## C2) AgeCl CANNOT CHANGE
         ## NOTE: This is because in the *new* approach AgeCl actually
         ##       represents COHORT CLASSES
   if ("AgeCl1" %in% vars) {
         Cmask[ Cmask$AgeCl1 != Cmask$AgeCl2, "Freq"] <- 1
    }

  # Expand the migration flows table using the count column 'Freq' in order to later
  # use ftable
  CCmask <- Cmask[rep(1:nrow(Cmask), Cmask$Freq), names(Cmask) != "Freq"]

  # Compute the flat contingency table
  CmaskT <- ftable(CCmask, row.vars = row.vars, col.vars = col.vars)

  # Now map cell values: 2 -> 0
  CmaskT[CmaskT == 2] <- 0

  # If needed print
  if (print) print(CmaskT, zero.print = ".")

  # Return output object
  return(CmaskT)
}


`mkCheck` <- function(F, Cmask, verbose = TRUE) {
################################################
# CHECK that Migration Flows data are FEASIBLE #
################################################

  TEST <- F * Cmask
  UNF <- (TEST > 0)
  NUNF <- sum(UNF)
  if (NUNF > 0) {
      cat("\n##")
      cat("\n## WRONG migration flows: ", NUNF, "UNFEASIBLE cells detected!")
      cat("\n##\n")

      TEST.df <- as.data.frame(TEST)
      out <- TEST.df[TEST.df$Freq > 0, ]
      if (verbose) print(out)
    }
  else {
      cat("\n##")
      cat("\n## WELL FORMED migration flows: no unfeasible cells detected!")
      cat("\n##\n\n")
      out <- NULL
    }

  return(invisible(out))
}


`settle.AbroadP1` <- function(P1, setval = 0) {
#############################################################################
# If stocks and flows involve classification variable 'Region', then settle #
# the *initial* abroad population count (i.e. Region = 'E' in P1).          #
# By convention, this is done by setting this unknown value of P1 to a      #
# *conventional* value, controlled by argument 'setval'.                    #
#                                                                           #
# NOTE: By default 'setval' is set to 0. A viable alternative would be to   #
#       set it to a *very large* value, like 1E12.                          #
#                                                                           #
# NOTE: By properly settling also P2, the abroad region will enter the      #
#       equations seamlessly, but it WILL NOT IMPLY ANY BALANCING           #
#       DISCREPANCY IN ITSELF.                                              #
#############################################################################
  if ("Region" %in% names(P1)) {
     P1$N[P1$Region == 'E'] <- setval
    }
  return(P1)
}


`settle.AbroadP2` <- function(P2, N, M) {
#############################################################################
# If stocks and flows involve classification variable 'Region', then settle #
# the *final* abroad population count (i.e. Region = 'E' in P2).            #
# This is done by setting this unknown value of P2 to the computed value of #
# the r.h.s. of the DBE.                                                    #
#                                                                           #
# NOTE: If not explicitly specified, the abroad counts of Births, Deaths    #
#       and Natural Increase are automatically set to zero while reading    #
#       the input data (in Standard Format).                                #
#                                                                           #
# NOTE: This way the abroad region enters the equations seamlessly, but it  #
#       DOES NOT IMPLY ANY BALANCING DISCREPANCY IN ITSELF.                 #
#############################################################################
  if ("Region" %in% names(P2)) {
     P2$N[P2$Region == 'E'] <- P1$N[P2$Region == 'E'] + N$N[P2$Region == 'E'] + rowSums(M)[P2$Region == 'E']
    }
  return(P2)
}




#########################################################################
# S T A N D A R D    I N P U T    ---->    B A L A N C E R    I N P U T #
#########################################################################

`makeObj` <- function(P1, P2, B, D, N, F, M) {
######################################################
# Transform input STOCKS and FLOWS into objects that #
# can be fed to function Balance.                    #
######################################################
         # Build the list of objects to be balanced
         RP1 <- rbind(P1$N, deparse.level = 0)
         RP2 <- rbind(P2$N, deparse.level = 0)
         RB  <- rbind(B$N, deparse.level = 0)
         RD  <- rbind(D$N, deparse.level = 0)
         RN  <- rbind(N$N, deparse.level = 0)
         F  <- as.matrix(F); dimnames(F) <- NULL
         M  <- as.matrix(M); dimnames(M) <- NULL
         # Now transpose
         CP1 <- cbind(P1$N, deparse.level = 0)
         CP2 <- cbind(P2$N, deparse.level = 0)
         CB  <- cbind(B$N, deparse.level = 0)
         CD  <- cbind(D$N, deparse.level = 0)
         CN  <- cbind(N$N, deparse.level = 0)
         FT  <- t(F)
         MT  <- t(M)

         obj <-list(
                     RP1 = RP1,
                     RP2 = RP2,
                     RB  = RB,
                     RD  = RD,
                     RN  = RN,
                     F   = F,
                     M   = M,
                   # Now transpose
                     CP1 = CP1,
                     CP2 = CP2,
                     CB  = CB,
                     CD  = CD,
                     CN  = CN,
                     FT  = FT,
                     MT  = MT
                    )
         attr(obj, "schema") <- P1[, names(P1) != "N"]
     # Return
     return(obj)
    }


`makeObjV` <- function(obj, Cmask, varmodel = c("poisson", "poisson-skellam", "homoskedastic"), M.N.adj.zeros = FALSE,
                       vP2 = 1, vF = 1, vM = 1, vP1 = 0, vB = 0, vD = 0, vN = 0) {
##########################################################################################
# This functions builds variance matrices for the objects to be balanced.                #
#                                                                                        #
# varmodel = "poisson"         -> This option (the default one) assumes variances to be  #
#                                 proportional to input *absolute* values, as would be   #
#                                 the case for Poisson counts.                           #
#                                                                                        #
# varmodel = "poisson-skellam" -> As the "poisson" option above, *except for*            #
#                                  - Natural Increase N = B - D                          #
#                                  - Net Migrations   M = t(F) - F                       #
#                                 whose variances are assumed to be proportional to      #
#                                 the *sum* of the (positive) values entering the        #
#                                 differences, as would be the case for Skellam          #
#                                 variables.                                             #
#                                                                                        #
# varmodel = "homoskedastic"   -> All entries of a given object share the same variance  #
#                                 value (a constant). Variances of different objects     #
#                                 can differ.                                            #
#                                                                                        #
# NOTE: Argument M.N.adj.zeros is only relevant for *varmodel = "poisson-skellam"* and   #
#       only affects the way variances of N (Natural Increase) and M (Net migrations)    #
#       are modeled:                                                                     #
#       - M.N.adj.zeros = FALSE -> default option. Variances of 0 raw counts in M and N  #
#                                  are set to 0 even if this is not actually implied by  #
#                                  the Skellam distribution. These counts will *not* be  #
#                                  adjusted in the balancing process, thus behaving as   #
#                                  it would happen for varmodel = "poisson" (see below). #
#                                                                                        #
#       - M.N.adj.zeros = TRUE ->  variances of 0 raw counts in M and N are modeled as   #
#                                  truly Skellam, thus will generally be *non 0*.        #
#                                  Potentially, these counts will be adjusted in the     #
#                                  balancing process, at odds with what would happen for #
#                                  varmodel = "poisson" (see below).                     #
#                                                                                        #
# NOTE: In *MOST CASES* (i.e. regardless the specific 'varmodel') values whose raw       #
#       estimates are exactly 0 will have *0 variance* and thus *remain 0* also after    #
#       balancing. While this is implicit in the "poisson" varmodel, it can be (and, by  #
#       default, *is*) enforced also for the other varmodel options (see argument        #
#       'M.N.adj.zeros' above). This way:                                                #
#       - THE BALANCING PROCESS TREATS ALL RAW ZERO VALUES AS *STRUCTURAL ZEROS*: THEY   #
#         WILL NEVER BE CHANGED.                                                         #
#       - THE BALANCING PROCESS *DOES NOT INVENT NEW STOCKS AND FLOWS*, RATHER JUST      #
#         *ADJUSTS EXISTING STOCKS AND FLOWS*.                                           #
#                                                                                        #
# NOTE: The vObj parameters (e.g. vP2, vF, ...) are used to switch on (if vObj = 1) or   #
#       switch off (if vObj = 0) the variance of the corresponding object.               #
#       Consistently with the note above, IF THE VARIANCE OF AN OBJECT IS SWITCHED OFF,  #
#       THE BALANCING PROCESS WILL NOT MODIFY THE CORRESPONDING RAW ESTIMATES.           #
#                                                                                        #
##########################################################################################

  schema <- attr(obj, "schema")
  Cmask  <- as.matrix(Cmask); dimnames(Cmask) <- NULL
  varmodel <- match.arg(varmodel)

# Update obj to set variances
## P1
   if (varmodel == "poisson" || varmodel == "poisson-skellam") {
     obj$CP1 <- abs(obj$CP1) * vP1
    } else {
     obj$CP1 <- abs(sign(obj$CP1)) * vP1
    }
   if ("AgeCl" %in% names(schema)) {
     obj$CP1[schema$AgeCl == 'N', ] <- 0
    }
   if ("Region" %in% names(schema)) {
     obj$CP1[schema$Region == 'E', ] <- 0
    }
## B
   if (varmodel == "poisson" || varmodel == "poisson-skellam") {
     obj$CB <- abs(obj$CB) * vB
    } else {
     obj$CB <- abs(sign(obj$CB)) * vB
    }
   if ("AgeCl" %in% names(schema)) {
     obj$CB[schema$AgeCl != 'N', ] <- 0
    }
   if ("Region" %in% names(schema)) {
     obj$CB[schema$Region == 'E', ] <- 0
    }
## D
   if (varmodel == "poisson" || varmodel == "poisson-skellam") {
     obj$CD <- abs(obj$CD) * vD
    } else {
     obj$CD <- abs(sign(obj$CD)) * vD
    }
   if ("Region" %in% names(schema)) {
     obj$CD[schema$Region == 'E', ] <- 0
    }
## N
   if (varmodel == "poisson") {
     # Traditional assumption: treat abs(N) as a poisson
     obj$CN <- abs(obj$CN) * vN
    } else  {
         if (varmodel == "poisson-skellam") {
         #  ----------------------------------------------------------------------------------------------------------  #
         #  # NEW assumption: treat N = B - D as a difference of two poissons (following the SKELLAM distribution)!
             ## IF raw zeros are not allowed to be adjusted, SET THEIR VARIANCE TO ZERO
             if (!M.N.adj.zeros) N_zeros <- !(obj$CN != 0)

             ## SET SKELLAM VARIANCE
             obj$CN <- ( obj$CB + obj$CD ) * vN

             ## IF raw zeros are not allowed to be adjusted, SET THEIR VARIANCE TO ZERO
             if (!M.N.adj.zeros) obj$CN[N_zeros] <- 0
         #  ----------------------------------------------------------------------------------------------------------  #
            } else {
             obj$CN <- abs(sign(obj$CN)) * vN
            }
        }
   if ("Region" %in% names(schema)) {
     obj$CN[schema$Region == 'E', ] <- 0
    }
## F
   if (varmodel == "poisson" || varmodel == "poisson-skellam") {
     obj$F <- abs(obj$F) * vF
    } else {
     obj$F <- abs(sign(obj$F)) * vF
    }
   obj$F <- obj$F * (1 - Cmask)
   diag(obj$F) <- 0
## M
   if (varmodel == "poisson") {
     # Traditional assumption: treat abs(M) as a poisson
     obj$M <- abs(obj$M) * vM
    } else  {
         if (varmodel == "poisson-skellam") {
         #  ----------------------------------------------------------------------------------------------------------  #
         #  # NEW assumption: treat M = t(F) - F as a difference of two poissons (following the SKELLAM distribution)!
             ## IF raw zeros are not allowed to be adjusted, SET THEIR VARIANCE TO ZERO
             if (!M.N.adj.zeros) M_zeros <- !(obj$M != 0)

             ## SET SKELLAM VARIANCE
             obj$M <- ( t(obj$F) + obj$F ) * vM

             ## IF raw zeros are not allowed to be adjusted, SET THEIR VARIANCE TO ZERO
             if (!M.N.adj.zeros) obj$M[M_zeros] <- 0
         #  ----------------------------------------------------------------------------------------------------------  #
            } else {
             obj$M <- abs(sign(obj$M)) * vM
            }
        }
   diag(obj$M) <- 0
## P2
   if (varmodel == "poisson" || varmodel == "poisson-skellam") {
     obj$CP2 <- abs(obj$CP2) * vP2
    } else {
     obj$CP2 <- abs(sign(obj$CP2)) * vP2
    }
#   # DON'T DO THAT WITH COHORTS!!!!!!!!
#   if ("AgeCl" %in% names(schema)) {
#     obj$CP2[schema$AgeCl == 'N', ] <- 0
#    }

# Now transpose
   obj$RP1 <- t(obj$CP1)
   obj$RB  <- t(obj$CB)
   obj$RD  <- t(obj$CD)
   obj$RN  <- t(obj$CN)
   obj$FT  <- t(obj$F)
   obj$MT  <- t(obj$M)
   obj$RP2 <- t(obj$CP2)

  attr(obj, "schema") <- schema
  # Return
  return(obj)
}




##################################
# S T A N D A R D    O U T P U T #
##################################

`write.ALL` <- function(BalObj, STOCK, FLOWS,
                        path = if (interactive()) choose.dir(getwd(), caption = "SELECT OUTPUT FOLDER") else NA,
                        file.ext = "csv", sep = ";", dec = '.', ...){
##############################################################################
# Convert balanced objects to STANDARD OUTPUT FORMAT and write them to       #
# external files:                                                            #
# - POPULATION COUNTS: QP1, QP2                                              #
# - BIRTHS:            QB                                                    #
# - DEATHS:            QD                                                    #
# - NATURAL INCREASE:  QN                                                    #
# - MIGRATION FLOWS:   QF                                                    #
# - NET MIGRATIONS:    QM                                                    #
#                                                                            #
# NOTE: BalObj must be a return object of function Balance.                  #
#                                                                            #
# NOTE: STOCK must be an already prepared POPULATION COUNTS object.          #
#       FLOWS must be an already prepared MIGRATION FLOWS object.            #
#       Both serve the purpose to guarantee stocks-flows FORMAT consistency  #
#       (thus only their *structure* matters, not their *counts*).           #
#                                                                            #
# NOTE: If 'path' is NA, no external files are written.                      #
#                                                                            #
# REMARK: The purpose of this function is reporting. For parsimony, output   #
#         objects are built keeping only meaningful modalities.              #
#         Therefore output objects will *not* be congruent with their input  #
#         counterparts (as returned by read.P(), read.BDN() and read.F()).   #
##############################################################################

  # Get balanced objects names
  BO.names <- names(BalObj)

  # Is BalObj balanced or not?
  # NOTE: The test is FALSE even if BalObj *is balanced* WHEN it HAS BEEN SET
  #       against the constraints via function *set.BalObj*. This causes no
  #       troubles, as names.out are correctly built through the if clauses
  #       below
  is.balanced <- identical(unique(substr(BO.names, 1, 1)), "Q")

  # POPULATION COUNTS, BIRTHS, DEATHS, NATURAL INCREASE
  PN.names <- BO.names[(endsWith(BO.names, "CP1") |  endsWith(BO.names, "CP2") | endsWith(BO.names, "CB")  | endsWith(BO.names, "CD") | endsWith(BO.names, "CN"))]
  PN.names.out <- gsub("C", "", PN.names)
  if (!is.balanced) {
     PN.names.out <- paste("Q", PN.names.out, sep = "")
    }

  # Initialize output objects list
  OUTlist <- NULL
  
  for (obj.name in PN.names){
       # Add the right "schema" to the balanced output
       obj <- STOCK
       # Add the column of balanced figures
       obj$N <- as.numeric(BalObj[[obj.name]])

       # Restrict to the meaningful modalities:
       # POPULATION COUNTS,
       # BIRTHS, DEATHS, NATURAL INCREASE -> exclude 'abroad'  level 'E' of Region
       # POPULATION COUNTS                -> exclude 'newborn' level 'N' of AgeCl       # DON'T DO THAT WITH COHORTS!!!!!!!!!!!!
       # BIRTHS                           -> exclude all AgeCl levels but 'N' of AgeCl
       if ("Region" %in% names(obj)) {
         obj <- obj[obj$Region != "E", ]
        }
#       if ("AgeCl" %in% names(obj) && ( endsWith(obj.name, "CP1") || endsWith(obj.name, "CP2") ) ) {
#         obj <- obj[obj$AgeCl != "N", ]
#        }
       if ("AgeCl" %in% names(obj) && endsWith(obj.name, "CB") ) {
         obj <- obj[obj$AgeCl == "N", ]
        }

       rownames(obj) <- NULL

       # Add object to output list
       ol <- list(obj)
       names(ol) <- gsub("Q", "", PN.names.out[PN.names == obj.name])
       OUTlist <- c(OUTlist, ol)

       # Write to external file
       if (!is.na(path)) {
         write.table(x = obj,
                     file = file.path(path, paste(PN.names.out[PN.names == obj.name], file.ext, sep=".")),
                     quote = FALSE,
                     sep = sep,
                     dec = dec,
                     row.names = FALSE,
                     col.names = TRUE,
                     ...)
        }
    }

  # MIGRATION FLOWS, NET MIGRATIONS
  FM.names <- BO.names[(endsWith(BO.names, "F") |  endsWith(BO.names, "M"))]
  FM.names.out <- FM.names
  if (!is.balanced) {
     FM.names.out <- paste("Q", FM.names.out, sep = "")
    }

  for (obj.name in FM.names){
       # Generate the right cols and rows metadata
       obj <- BalObj[[obj.name]] * (1 - FLOWS + FLOWS)
       # Convert to data.frame format
       obj <- as.data.frame(obj)
       # Restrict to Freq != 0 entries when reporting (save space and downstream analysis burden)
       obj <- obj[obj$Freq != 0, ]

       names(obj)[names(obj) == "Freq"] <- "N"

       # Order according to factors
       obj.factors <- names(obj)[names(obj) != "N"]
       obj <- obj[do.call(order, obj[, obj.factors]), ]

       # NO NEED to restrict to the meaningful modalities: non-zero counts
       # can only arise from meaningful flows!
       rownames(obj) <- NULL

       # Add object to output list
       ol <- list(obj)
       names(ol) <- gsub("Q", "", FM.names.out[FM.names == obj.name])
       OUTlist <- c(OUTlist, ol)

       # Write to external file
       if (!is.na(path)) {
         write.table(x = obj,
                     file = file.path(path, paste(obj.name, file.ext, sep=".")),
                     quote = FALSE,
                     sep = sep,
                     dec = dec,
                     row.names = FALSE,
                     col.names = TRUE,
                     ...)
        }
    }

  return(invisible(OUTlist))
}


`compare.ALL` <- function(BalObj, Obj, STOCK, FLOWS,
                          path = if (interactive()) choose.dir(getwd(), caption = "SELECT OUTPUT FOLDER") else NA,
                          file.ext = "csv", sep = ";", dec = '.', ...){
########################################################################
# Convert initial and balanced objects to STANDARD OUTPUT FORMAT,      #
# compare homologous values and write to external files.               #
#                                                                      #
# REMARK: The purpose of this function is reporting. For parsimony,    #
#         output objects are built keeping only meaningful modalities. #
#         Therefore output objects will *not* be congruent with their  #
#         input counterparts (as returned by read.P(), read.BDN() and  #
#         read.F()).                                                   #
########################################################################

  # Is BalObj balanced or not?
  BalObj.balanced <- identical(unique(substr(names(BalObj), 1, 1)), "Q")
  if (!BalObj.balanced) stop("BalObj must contain *balanced* objects!")

  # Is Obj balanced or not?
  Obj.balanced <- identical(unique(substr(names(Obj), 1, 1)), "Q")
  if (Obj.balanced) stop("Obj must contain *unbalanced* objects!")

  # Original
  O <- write.ALL(Obj, STOCK, FLOWS, path = NA, file.ext, sep, dec, ...)
  # Sort according to a conventional COMMON order
  O <- O[sort(names(O))]

  # Balanced
  Q <- write.ALL(BalObj, STOCK, FLOWS, path = NA, file.ext, sep, dec, ...)
  # Sort according to a conventional COMMON order
  Q <- Q[sort(names(O))]

  # Comparison: O vs. Q
  ### OLD WAY: Faster, but would *not work* when: ( varmodel = "poisson-skellam" ) AND ( M.N.adj.zeros = TRUE )
  # COMP <- mapply(FUN = function(o, q) cbind(o, Nbal = q$N), O, Q)

  ### NEW WAY: Slower, but *general*
  O.Q.merge <- function(o, q) {
     names(q)[names(q) == "N"] <- "Nbal"
     O.Q <- merge(o, q, sort = TRUE, all = TRUE)
     O.Q$N[is.na(O.Q$N)] <- 0
     O.Q$Nbal[is.na(O.Q$Nbal)] <- 0
     # Order according to factors
     O.Q.factors <- names(O.Q)[!(names(O.Q) %in% c("N", "Nbal"))]
     O.Q <- O.Q[do.call(order, O.Q[, O.Q.factors]), ]
     rownames(O.Q) <- NULL
     return(O.Q)
    }
  COMP <- mapply(FUN = O.Q.merge, O, Q)
  names(COMP) <- paste("Delta", names(COMP), sep = "")


  # Write to external file
  if (!is.na(path)) {
     for (obj.name in names(COMP)){
         obj <- COMP[[obj.name]]
         write.table(x = obj,
                     file = file.path(path, paste(obj.name, file.ext, sep=".")),
                     quote = FALSE,
                     sep = sep,
                     dec = dec,
                     row.names = FALSE,
                     col.names = TRUE,
                     ...)
        }
    }

  return(invisible(COMP))
}




########################################
# O U T P U T    D I A G N O S T I C S #
########################################


`getMargins` <- function(Obj, P1, P2, B, D, N, F, M, margin) {
###############################################################
# This function marginalizes stocks and flows values in Obj   #
# with respecto to classification variable 'margin'.          #
#                                                             #
# REMARK: Aggregation on margins here *includes convenience   #
#         values that are absent* from outputs of functions   #
#         write.ALL (Q) and compare.ALL (Delta). These can    #
#         sometimes be *NON MEANINGFUL*, e.g. the value of P2 #
#         for Region == "E" (see function settle.AbroadP2).   # 
#                                                             #
# NOTE: DBE errors at aggregated level are also reported.     #
###############################################################

  # Get balanced objects names
  O.names <- names(Obj)

  # Is Obj balanced or not?
  # NOTE: The test is FALSE even if Obj *is balanced* WHEN it HAS BEEN SET
  #       against the constraints via function *set.BalObj*. This causes no
  #       troubles, as names.out are correctly built through the if clauses
  #       below
  is.balanced <- identical(unique(substr(O.names, 1, 1)), "Q")

  if (!is.balanced) {
     # For convenience, Obj names will *always* start with 'Q' *even* when Obj is *not balanced* 
     names(Obj) <- paste("Q", O.names, sep = "")
    }

  # Check margin variable is just one and exists
  if (length(margin) > 1) stop("Please specify just one margin variable!")
  if (!(margin %in% names(P1)[names(P1) != "N"])) stop("Illicit margin variable!")

  # Build marginalization formula
  margin.formula <- as.formula(paste("N ~", margin), env = .GlobalEnv)

  # Build and aggregate DBE terms: P1, P2, B, D, N, F, M
  ## P1
  QP1 <- P1
  QP1$N <- as.numeric(Obj$QCP1)
  P1 <- aggregate(margin.formula, data = QP1, FUN = sum)

  ## P2
  QP2 <- P2
  QP2$N <- as.numeric(Obj$QCP2)
  P2 <- aggregate(margin.formula, data = QP2, FUN = sum)

  ## B
  QB <- B
  QB$N <- as.numeric(Obj$QCB)
  B <- aggregate(margin.formula, data = QB, FUN = sum)

  ## D
  QD <- D
  QD$N <- as.numeric(Obj$QCD)
  D <- aggregate(margin.formula, data = QD, FUN = sum)

  ## N
  QN <- N
  QN$N <- as.numeric(Obj$QCN)
  N <- aggregate(margin.formula, data = QN, FUN = sum)

  # Build row and col vars for flat contingency table marginalization
  row.margin <- paste(margin, "1", sep = "")
  col.margin <- paste(margin, "2", sep = "")
  
  ## F
  QF <- Obj$QF * (1 - F + F)
  F <- ftable(QF, row.vars = row.margin, col.vars = col.margin)
  F <- zapsmall(F)

  ## M
  QM <- Obj$QM * (1 - M + M)
  M <- ftable(QM, row.vars = row.margin, col.vars = col.margin)
  M <- zapsmall(M)

  # Compute aggregated DBE errors
  Err <- P2
  names(Err)[names(Err) == "N"] <- "P2"
  Err$P1 <- P1$N
  Err$N <- N$N
  Err$rowSumsM <- rowSums(M)
  Err$ERR <- Err$P2 - Err$P1 - Err$N - Err$rowSumsM
  Err[, -1] <- zapsmall(as.matrix(Err[, -1]))

  # Build output list
  Margins <- list(P1 = P1, P2 = P2, B = B, D = D, N = N, F = F, M = M)

  attr(Margins, "DBEerrors") <- Err
  attr(Margins, "margin") <- margin
  attr(Margins, "is.balanced") <- is.balanced

  return(Margins)
}


`getDBEerrors` <- function(Obj, P1, P2, N, M, tol = 1E-7, verbose = TRUE) {
###############################################################################
# This function checks for residual discrepancies in atomic DBEs after the    #
# balancing process and returns relevant diagnostic information.              #
#                                                                             #
# NOTE: If Obj is *not balanced* atomic DBE discrepancies refer to raw data.  #
#       In this case, a large number of discrepancies is expected and         #
#       returned diagnostics will likely be heavy.                            #
#                                                                             #
# NOTE: The number of non satisfied atomic DBEs reported by this function is  #
#       the *actual* one, i.e. *for raw data* HALF the number of broken micro #
#       constraints detected by workhorse function Balance(). Indeed,         #
#       Balance() - for numerical stability and computational efficiency      #
#       considerations - uses TWICE the actual DBEs.                          #
###############################################################################

  # Get balanced objects names
  O.names <- names(Obj)

  # Is Obj balanced or not?
  # NOTE: The test is FALSE even if Obj *is balanced* WHEN it HAS BEEN SET
  #       against the constraints via function *set.BalObj*. This causes no
  #       troubles, as names.out are correctly built through the if clauses
  #       below
  is.balanced <- identical(unique(substr(O.names, 1, 1)), "Q")

  if (!is.balanced) {
     # For convenience, Obj names will *always* start with 'Q' *even* when Obj is *not balanced* 
     names(Obj) <- paste("Q", O.names, sep = "")
    }

  # Build DBE terms: P1, P2, N, M
  ## P1
  QP1 <- P1
  QP1$N <- as.numeric(Obj$QCP1)
  ## P2
  QP2 <- P2
  QP2$N <- as.numeric(Obj$QCP2)
  ## N
  QN <- N
  QN$N <- as.numeric(Obj$QCN)
  ## M
  QM <- Obj$QM * (1 - F + F)

  #######################
  # Atomic level checks #
  #######################

  # Number of atomic DBEs
  nDBEs <- NROW(P1)

  # Compute DBE discrepancies at micro level
  rowSumsM <- rowSums(QM)  # Just to reuse it later
  E <- QP2$N - QP1$N - QN$N - rowSumsM

  # Compute the Total Absolute Error (numerically zero as it SHOULD BE?)
  TotAbsError <- sum(abs(E))

  # Compute the Total Error,  i.e. the population level DBE discrepancy (numerically zero as it SHOULD BE?)
  TotError <- sum(E)

  # Discrepancies at tolerance level tol
  isErr <- (abs(E) > tol)

  if (any(isErr)) {

     # Build output dataframe
     Err <- QP2[isErr, ]
     names(Err)[names(Err) == "N"] <- "P2"
     Err$P1 <- QP1$N[isErr]
     Err$N <- QN$N[isErr]
     Err$rowSumsM <- rowSumsM[isErr]
     Err$ERR <- E[isErr]
     # Number of broken atomic DBEs at tolerance level tol
     nDBEs.KO <- NROW(Err)
     attr(Err, "nDBEs") <- nDBEs
     attr(Err, "tol") <- tol
     attr(Err, "nDBEs.KO") <- nDBEs.KO
     attr(Err, "TotAbsError") <- TotAbsError
     attr(Err, "TotError") <- TotError

     if (verbose) {
         # Print some diagnostic messages
         cat("- Atomic DBEs not satisfied (at tolerance level ", tol, "): ", nDBEs.KO, " (out of ", nDBEs, ")\n", sep = "")
         cat("- Total Absolute Error: ", TotAbsError, "\n", sep = "")
         cat("- Total Error: ", TotError, "\n\n", sep = "")
         cat("- Errors dataframe (up to first 30 rows):", "\n\n", sep = "")
         print(head(Err, 30))
        }

    } else {

     # Build output string
     Err <- paste("All atomic DBEs satisfied (at tolerance level ", tol, ")", sep = "")
     # Number of broken atomic DBEs at tolerance level tol
     nDBEs.KO <- 0
     attr(Err, "nDBEs") <- nDBEs
     attr(Err, "tol") <- tol
     attr(Err, "nDBEs.KO") <- nDBEs.KO
     attr(Err, "TotAbsError") <- TotAbsError
     attr(Err, "TotError") <- TotError

     if (verbose) {
         # Print some diagnostic messages
         cat("- ", Err, "\n", sep = "")
         cat("- Total Absolute Error: ", TotAbsError, "\n", sep = "")
         cat("- Total Error: ", TotError, "\n\n", sep = "")
        }
    }

  # Return output object invisibly
  class(Err) <- c("DBEerrors", class(Err))
  return(invisible(Err))
}

`print.DBEerrors` <- function(x, ...) {
  nDBEs <- attr(x, "nDBEs")
  tol <- attr(x, "tol")
  nDBEs.KO <- attr(x, "nDBEs.KO")
  TotAbsError <- attr(x, "TotAbsError")
  TotError <- attr(x, "TotError")
  # Print some diagnostic messages

  if (!is.character(x)) {
     # SOME broken atomic DBEs at tolerance level tol
     cat("- Atomic DBEs not satisfied (at tolerance level ", tol, "): ", nDBEs.KO, " (out of ", nDBEs, ")\n", sep = "")
     cat("- Total Absolute Error: ", TotAbsError, "\n", sep = "")
     cat("- Total Error: ", TotError, "\n\n", sep = "")
     cat("- Errors dataframe (up to first 30 rows):", "\n\n", sep = "")
     print.data.frame(head(x, 30))
    } else {
     # NO broken atomic DBEs at tolerance level tol
     cat("- ", x, "\n", sep = "")
     cat("- Total Absolute Error: ", TotAbsError, "\n", sep = "")
     cat("- Total Error: ", TotError, "\n\n", sep = "")
    }

  return(invisible(x))
}


`checkNegativeEstimates` <- function(Delta, test = c("bal", "raw")) {
########################################################################
# This functions checks for *negative* balanced estimates of counts    #
# and returns relevant diagnostic information.                         #
#                                                                      #
# NOTE: Delta must be the output of comparison function compare.ALL(). #
#                                                                      #
# NOTE: If test = "raw" the check is performed on raw data. In this    #
#       case, *zero* negative counts are expected.                     #
########################################################################

  # Get balanced objects names
  test <- match.arg(test)

  CountVar <- switch(test, bal = "Nbal", raw = "N")

  # Check counts that should be positive: P1, P2, B, D, F
  ## P1
  P1 <- Delta$DeltaP1[Delta$DeltaP1[, CountVar] < 0, ]
  ## P2
  P2 <- Delta$DeltaP2[Delta$DeltaP2[, CountVar] < 0, ]
  ## B
  B <- Delta$DeltaB[Delta$DeltaB[, CountVar] < 0, ]
  ## D
  D <- Delta$DeltaD[Delta$DeltaD[, CountVar] < 0, ]
  ## F
  F <- Delta$DeltaF[Delta$DeltaF[, CountVar] < 0, ]

  # Build output list
  Neg <- list(P1 = P1, P2 = P2, B = B, D = D, F = F)

  return(Neg)
}


`MicroMonitor` <- function(Delta,
                           i = c(Sex1 = NA, Nation1 = NA, AgeCl1 = NA, Region1 = NA),
                           j = c(Sex2 = NA, Nation2 = NA, AgeCl2 = NA, Region2 = NA)
                        ) {
###################################################################################
# Given input profiles i and (optionally) j, the objects entering the micro-level #
# DBE:                                                                            #
#                                                                                 #
# DBE_i -> P2_i = P1_1 + (B_i - D_i) + SUM_j(F_ji - F_ij) =                       #
#               = P1_1 +     N_i     + SUM_j(    M_ij   )                         #
#                                                                                 #
# are returned.                                                                   #
#                                                                                 #
# NOTE: Delta must be the output of comparison function compare.ALL()             #
#                                                                                 #
# NOTE: if j is NOT passed, the 2-index objects above are returned for *all the   #
#       relevant values of j*. Otherwise, the result is restricted *only to the   #
#       specified j*.                                                             # 
###################################################################################
  if (is.null(names(i))) names(i) <- c("Sex1", "Nation1", "AgeCl1", "Region1")

  # P1 population counts
  P1_ij <- with(Delta$DeltaP1, Delta$DeltaP1[Sex == i["Sex1"] & Nation == i["Nation1"] & AgeCl == i["AgeCl1"] & Region == i["Region1"], ] )

  # P2 population counts
  P2_ij <- with(Delta$DeltaP2, Delta$DeltaP2[Sex == i["Sex1"] & Nation == i["Nation1"] & AgeCl == i["AgeCl1"] & Region == i["Region1"], ] )

  # Births: expect them to be NULL
  B_ij <- with(Delta$DeltaB, Delta$DeltaB[Sex == i["Sex1"] & Nation == i["Nation1"] & AgeCl == i["AgeCl1"] & Region == i["Region1"], ] )

  # Deaths
  D_ij <- with(Delta$DeltaD, Delta$DeltaD[Sex == i["Sex1"] & Nation == i["Nation1"] & AgeCl == i["AgeCl1"] & Region == i["Region1"], ] )

  # Natural Increase
  N_ij <- with(Delta$DeltaN, Delta$DeltaN[Sex == i["Sex1"] & Nation == i["Nation1"] & AgeCl == i["AgeCl1"] & Region == i["Region1"], ] )

  # Migrations
  if (missing(j)) {
     # Inflows to i from j
     F_ji <- with(Delta$DeltaF, Delta$DeltaF[Sex2 == i["Sex1"] & Nation2 == i["Nation1"] & AgeCl2 == i["AgeCl1"] & Region2 == i["Region1"], ] )
     # Outflows from i to j
     F_ij <- with(Delta$DeltaF, Delta$DeltaF[Sex1 == i["Sex1"] & Nation1 == i["Nation1"] & AgeCl1 == i["AgeCl1"] & Region1 == i["Region1"], ] )
     # Append Inflows to i and outflows from i
     F_ij <- rbind(F_ji, F_ij)
     # Net Migrations
     M_ij <- with(Delta$DeltaM, Delta$DeltaM[Sex1 == i["Sex1"] & Nation1 == i["Nation1"] & AgeCl1 == i["AgeCl1"] & Region1 == i["Region1"], ] )
    } else {
     if (is.null(names(j))) names(j) <- c("Sex2", "Nation2", "AgeCl2", "Region2")
     # Inflows to i from j
     F_ji <- with(Delta$DeltaF, Delta$DeltaF[Sex1 == j["Sex2"] & Nation1 == j["Nation2"] & AgeCl1 == j["AgeCl2"] & Region1 == j["Region2"] &
                                             Sex2 == i["Sex1"] & Nation2 == i["Nation1"] & AgeCl2 == i["AgeCl1"] & Region2 == i["Region1"], ] )
     # Outflows from i to j
     F_ij <- with(Delta$DeltaF, Delta$DeltaF[Sex1 == i["Sex1"] & Nation1 == i["Nation1"] & AgeCl1 == i["AgeCl1"] & Region1 == i["Region1"] &
                                             Sex2 == j["Sex2"] & Nation2 == j["Nation2"] & AgeCl2 == j["AgeCl2"] & Region2 == j["Region2"], ] )
     # Append Inflows to i and outflows from i
     F_ij <- rbind(F_ji, F_ij)
     # Net Migrations
     M_ij <- with(Delta$DeltaM, Delta$DeltaM[Sex1 == i["Sex1"] & Nation1 == i["Nation1"] & AgeCl1 == i["AgeCl1"] & Region1 == i["Region1"] &
                                     Sex2 == j["Sex2"] & Nation2 == j["Nation2"] & AgeCl2 == j["AgeCl2"] & Region2 == j["Region2"], ] )
    }

  out <- list(P1 = P1_ij, P2 = P2_ij, B = B_ij, D = D_ij, N = N_ij, F = F_ij, M = M_ij)

  return(out)
}

