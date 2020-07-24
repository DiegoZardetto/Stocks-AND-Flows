#############################################################################
# A set of functions for solving National Accounts balancing problems.      #
#                                                                           #
# SW name:     BalanceR (could become an R package, eventually)             #
# Author:      Diego Zardetto                                               #
# SW stage:    Pre-alpha                                                    #
# Istat ref:   D08-182/DGEN-19/11/2013                                      #
# Credits:     This is an independent R implementation of V. Nicolardi's    #
#              original Speakeasy procedure.                                #
# Disclaimer:  This software is currently being developed for research      #
#              purposes, thus you could only use it at your own risk.       #
# Copying:     Except as otherwise stated, the code is copyright 2014       #
#              Diego Zardetto ©.                                            #
#############################################################################


`read.BalConstr` <- function(file = if (interactive()) choose.files(multi = FALSE, caption = "SELECT CONSTRAINTS FILE"),
                             skip.decl = TRUE, verbose = TRUE, time = FALSE, ...){
#########################################################################
# Read an account system expressed in *eulero* sintax from a text file. #
# NOTE: Returns a list of data structures describing the system.        #
# NOTE: Argument skip.decl enables/disables skipping declarations (i.e. #
#       fake constraints which are a relic of Nicolardi's syntax).      #
#       Since declarations are useless, they could be directly avoided  #
#       by the user when writing the constraint file. Should this be    #
#       not the case, they can be skipped by read.BalConstr. In either  #
#       case, balancing will be faster and provide identical results.   #
# NOTE: Argument 'time' and related elaborations are here only          #
#       temporarily and will be removed when Speakeasy Task Force is    #
#       over.                                                           #
#########################################################################

  # Read single constraints, one per line
  v.constr <- readLines(con = file, warn = FALSE)

  # If requested start timing
  if ( Time <- (time && !missing(file)) ) Tstart <- proc.time()

  getTags <- function(constr){
     # Split each constraint by "white space", i.e. one *or more* spaces
     tags <- unlist(strsplit(constr, "[[:space:]]{1,}"))
     # Remove empty strings (i.e. "") arising from initial spaces (if any)
     tags[nchar(tags) > 0]
    }
  # Constraints list: each component stores the tags of the corresponding constraint
  l.constr <- lapply(v.constr, getTags)
  n.constr <- length(l.constr)
  n.tag.constr <- sapply(l.constr, length)
  max.tag <- max(n.tag.constr)
  # Initialize the (character) constraint matrix
  constr.mat <- matrix(as.character(NA), n.constr, max.tag)
  # Fill the (character) constraint matrix
  sapply(1:n.constr, function(i) constr.mat[i, 1:n.tag.constr[i]] <<- l.constr[[i]])
  # Synctactically correct blocks involve 3 tags: (i) sign, (ii) operator, (iii) object -> this should be checked
  eq.blocks <- max.tag / 3
  sgn.col <- seq(1, by = 3, length.out = eq.blocks)
  op.col  <- seq(2, by = 3, length.out = eq.blocks)
  obj.col <- seq(3, by = 3, length.out = eq.blocks)
  # (i) Signs matrix (numeric, eventually)
  sgn.mat <- constr.mat[, sgn.col, drop=FALSE]
    sgn.mat <- matrix(c(-1, 1)[match(sgn.mat, c("-", "+"))], nrow(sgn.mat), ncol(sgn.mat))
    # If one wants to map NA signs (arising from missing tags) to 0
    sgn.mat[is.na(sgn.mat)] <- 0
  # (ii) Operators matrix (character + integer re-code)
  op.mat <- constr.mat[, op.col, drop=FALSE]
    # Initialize integer matrix (to be)
    op.code.mat <- op.mat
    # Map NAs to 0 (eventually)
    op.code.mat[is.na(op.code.mat)] <- 0
    # Code operators into integers
    op.code <- c("SR" = 1, "SC" = 2, "SM" = 3, "MM" = 4, "VC" = 4, "VR" = 4)
    # Recode op.mat
    sapply(names(op.code), function(op) op.code.mat[op.code.mat==op] <<- op.code[op])
    mode(op.code.mat) <- "numeric"
  # (iii) Object matrix (character)
  obj.mat <- constr.mat[, obj.col, drop=FALSE]

  # If required, identify declarations (fake constraints, inherited from
  # Nicolardi's syntax) and drop them (18/06/2014)
  if (skip.decl) {
      # A constraint is deemed to be a declaration when any obj appears in it
      # more than once (maybe some more robust rule could be devised...) 
      is.true.constr <- apply(obj.mat, 1, function(obj) !any(duplicated(obj, incomparables = NA)))
      # Restrict data structures to true constraints
      constr.mat <- constr.mat[is.true.constr, , drop = FALSE]
      sgn.mat <- sgn.mat[is.true.constr, , drop = FALSE]
      op.mat <- op.mat[is.true.constr, , drop = FALSE]
      op.code.mat <- op.code.mat[is.true.constr, , drop = FALSE]
      obj.mat <- obj.mat[is.true.constr, , drop = FALSE]
    }

  # Find the names of the object appearing in constraints: remove duplicates,
  # respect the order of appearance of constraints in the configuration file
  # and the one of objects within constraints
  # NOTE: In Nicolardi's code, since each object must appear at least into a
  #       "declaration" block, search is restricted to "non operational" blocks
  #       (op.code 4). THUS CONSTRAINTS MUST BE WRITTEN CAREFULLY, as objects
  #       appearing in the file but not tied to MM, VC or VR will generate
  #       INCONSISTENCIES!!!!!!!!!!!!

# DEBUG START: see NOTE below
#  obj.names <- unique(unlist(lapply(1:nrow(obj.mat), function(i.eq) obj.mat[i.eq, op.code.mat[i.eq, ] == 4])))

# NOTE: below a viable alternative permitting to neglect declarations.
# WARNING: NOT extensively tested yet!!!!!!!!! 
  obj.names <- sapply(1:nrow(obj.mat), function(i.eq) obj.mat[i.eq, ])
  obj.names <- unique(obj.names[!is.na(obj.names)])
# DEBUG END: see NOTE above

  n.obj <- length(obj.names)

  # Recode obj.mat by replacing each object with its position inside obj.names,
  # and mapping NAs into 0
  obj.pos.mat <- matrix(match(obj.mat, obj.names), nrow(obj.mat), ncol(obj.mat))
  obj.pos.mat[is.na(obj.mat)] <- 0

  # For each object, identify actual constraints (excluding daclarations) involving
  # it and the sign with which it enters
    # First find all eqs in which obj appears and where it appears in them
    constr.l <- lapply(obj.names, function(name) {
                                           m <- which(obj.mat==name, arr.ind=TRUE)
                                           m[order(m[, "row"]), , drop=FALSE]})
    # Get constraints positions list
    constr.pos.l <- lapply(constr.l, function(m) m[, "row"])

    # Then remove those eqs in which obj appears more than once (indeed, this
    # can happen in declarations only).
    # NOTE: there should be at least a TRUE -> should check and signal if there are objs who are not actually involved in any REAL constr...
    true.constr.pos.l <- lapply(constr.pos.l, function(eqs) !(eqs %in% eqs[duplicated(eqs)]))

    # Restrict constraints positions list to true constraints...
    constr.pos.l <- lapply(1:n.obj, function(i) constr.pos.l[[i]][true.constr.pos.l[[i]]])

    #... and do the same with the signs
    true.constr.l <- lapply(1:n.obj, function(i) constr.l[[i]][true.constr.pos.l[[i]], , drop=FALSE])
    constr.sgn.l  <- lapply(1:n.obj, function(i) {
                                              m <- true.constr.l[[i]]
                                              sapply(1:nrow(m), function(j.row) sgn.mat[m[j.row, "row"], m[j.row, "col"]]) })

    # Lastly put all into a matrix (NOTE: here number of columns is DIFFERENT from
    # [pointeq] and [signseq] of Francesco's code, which is flawed (because, e.g.,
    # an object *can* appear into more true eqs than macro blocks)
    n.constr.obj <- sapply(constr.pos.l, length)    
    max.constr.obj <- max(n.constr.obj)
    constr.pos.obj <- sgn.pos.obj <- matrix(0, n.obj, max.constr.obj)
    sapply(1:n.obj, function(i) constr.pos.obj[i, 1:n.constr.obj[i]] <<- constr.pos.l[[i]])
    sapply(1:n.obj, function(i) sgn.pos.obj[i, 1:n.constr.obj[i]] <<- constr.sgn.l[[i]])

    # Store into a vector the positions of true constraints
    true.constr.i <- sort(unique(as.numeric(constr.pos.obj[constr.pos.obj>0])))
    n.true.constr <- length(true.constr.i)

    # Text of actual constraints (no declarations)
	true.v.constr <- if (skip.decl) v.constr[is.true.constr] else v.constr[true.constr.i]

  # If requested stop timing
  if (Time) Tstop <- proc.time()

  # Build output list (must consider the opportunity of defining a class...)
  BalConstr <- list(
                     v.constr = v.constr,                     # serves functions Balance() and check.ZeroVarConstr()
                     sgn.mat = sgn.mat,
                     op.mat = op.mat,                         # currently unused
                     op.code.mat = op.code.mat,
                     obj.mat = obj.mat,                       # currently unused
                     obj.names = obj.names,
                     obj.pos.mat = obj.pos.mat,
                     constr.pos.obj = constr.pos.obj,
                     sgn.pos.obj = sgn.pos.obj,
                     true.constr.i = true.constr.i,
                     true.v.constr = true.v.constr,           # serves functions Balance() and check.ZeroVarConstr()
                     skip.decl = skip.decl                    # serves functions Balance() and check.ZeroVarConstr()
                    )

  if (Time) attr(BalConstr, "time") <- (Tstop - Tstart)

# If needed, summarize the balancing system
if (verbose){
     cat("- Number of (macro) constraints in the balancing system: ", n.constr, "\n", sep="")
     cat("  of which:\n")
     cat("  - actual (macro) constraints: ", n.true.constr, "\n", sep="")
     cat("  - object declarations:        ", n.constr - n.true.constr, "\n", sep="")
     cat("- Number of objects involved in the balancing system:    ", n.obj, "\n\n", sep="")
    }

# Return output list
  BalConstr
}


`read.BalObj` <- function(path = if (interactive()) choose.dir(getwd(), caption = "SELECT INPUT FOLDER"),
                          file.ext = "csv", sep = ";", dec = '.', verbose = FALSE, ...){
###################################################################
# Import the objects (i.e. numeric matrices) to be balanced by    #
# reading *all* the files with a given extension from a specified #
# directory.                                                      #
# NOTE: Returns a list storing the objects, with names derived    #
#       from those of the corresponding input files.              #
# NOTE: Can be used to read both *objects* and *variances*,       #
#       provided the corresponding files are kept in *separate*   #
#       directories.                                              #
###################################################################

  # Initialize the list of objects to undergo balancing
  BalObj <- list()
  for (f in list.files(path, pattern = paste("\\.", file.ext, "$", sep=""))){ 
       obj.name <- substr(f, start = 1, stop = nchar(f) - 4)
       BalObj[[obj.name]] <- as.matrix(read.table(file.path(path, f),
                                                  header = FALSE,
                                                  sep = sep,
                                                  dec = dec,
                                                  colClasses = "numeric",
                                                  ...)
                                    )
       dimnames(BalObj[[obj.name]]) <- NULL
    }

  if (verbose) {
      cat("- Number of objects read: ", length(BalObj), "\n\n")
    }

  BalObj
}


`set.BalObj` <- function(BalObj, BalConstr){
######################################################################
# Set BalObj so that it can be *safely* used in conjunction with     #
# BalConstr.                                                         #
# Object names are checked against constraints. If check is passed,  #
# a re-order version of BalObj is returned, so that its elements are #
# positionally tied to BalConstr$obj.names.                          #
# NOTE: Check allows for *a single starting letter mismatch* (this   #
#       is a legacy choice, as CN would store balanced object with   #
#       the old name preceded by a "Q", and variances with a "V").   #
#       When such a case happens, the names of the objects stored in #
#       the re-ordered version of BalObj are copied from those in    #
#       BalConstr$obj.names while the original ones are stored in a  #
#       dedicated attribute.                                         #
######################################################################

  if (!is.null(attr(BalObj, "o.names"))){
     return(BalObj)
    }

  c.obj.names <- BalConstr$obj.names
  o.obj.names <- names(BalObj)
  o.len <- length(o.obj.names)
  o.ord <- 1:o.len
  names(o.ord) <- o.obj.names
  if (length(c.obj.names) > o.len){
      stop("Input objects are less than those involved in balancing constraints!\n\n")   # Could (should) devise a more readable msg
    }

  if (all(c.obj.names %in% o.obj.names)){
      BalObj <- BalObj[c.obj.names]
      o.ord <- o.ord[c.obj.names]
      attr(BalObj, "o.names") <- c.obj.names
      names(attr(BalObj, "o.names")) <- o.ord
    }
  else {
      firstchar <- substr(o.obj.names[1], 1, 1)
      co.obj.names <- paste(firstchar, c.obj.names, sep="")
      if (!all(co.obj.names %in% o.obj.names)){
          stop("Mismatch in object names!\n\n")       # Could (should) devise a more readable msg
        }
      BalObj <- BalObj[co.obj.names]
      names(BalObj) <- c.obj.names
      o.ord <- o.ord[co.obj.names]
      attr(BalObj, "o.names") <- co.obj.names
      names(attr(BalObj, "o.names")) <- o.ord
    }

BalObj
}

`reset.BalObj` <- function(BalObj){
##################################################################
# If object names have been changed by function set.BalObj, they #
# are now reset to their original values. The original ordering  #
# of objects is also reconstructed.                              #
# NOTE: The sequence reset-set is an identity, i.e. it must hold #
#       always TRUE the following:                               #
# identical(BalObj, reset.BalObj(set.BalObj(BalObj, BalConstr))) #
##################################################################
  if (!is.null(attr(BalObj, "o.names"))){
      names(BalObj) <- attr(BalObj, "o.names")
      BalObj <- BalObj[order(as.integer(names(attr(BalObj, "o.names"))))]
      attr(BalObj, "o.names") <- NULL
    }

  BalObj
}


`ConstrMisfit` <- function(BalConstr, BalObj, i.constr = NULL, must.set = TRUE){
########################################################################
# Given a set of object to be balanced, assess the degree of violation #
# of the (macro) constraints of the account system.                    #
########################################################################

if (must.set) BalObj <- set.BalObj(BalObj, BalConstr)

# Default: get residuals on *all* constraints
n.constr <- nrow(BalConstr$obj.pos.mat)
if (is.null(i.constr)) {
     i.constr <- 1:n.constr
    }
# Map operator codes into appropriate R functions
op.fun <- function(op.code, obj) switch(op.code, `1`=cbind(rowSums(obj)), `2`=rbind(colSums(obj)), `3`=matrix(sum(obj)), `4`=obj)
# NOTE: it seems one could substitute *all* the matricial structures of BalConstr
#       (which involve 0s & NAs) with lists: it would be easier and lighter...
# NOTE: index below could be a basis for diagnostics/debug
# i.d <<- 0
misfit <- lapply(i.constr, function(i) {
                             # i.d <<- i.d + 1
                             obj.pos <- (BalConstr$obj.pos.mat[i, ] > 0)
                             objs <- BalConstr$obj.pos.mat[i, obj.pos]
                             sgns <- BalConstr$sgn.mat[i, obj.pos]
                             op.codes <- BalConstr$op.code.mat[i, obj.pos]
                             n.objs <- length(objs)
                             res <- Reduce(`+`, lapply(1:n.objs, function(j.obj) sgns[j.obj] *
                                                                          op.fun(op.codes[j.obj], BalObj[[objs[j.obj]]])))
                            })
misfit
}

`ConstrCumVar` <- function(BalConstr, BalObjV, i.constr = NULL, must.set = TRUE){
########################################################################
# Given the variances of a set of object to be balanced, cumulate the  #
# variances along the blocks of the (macro) constraints of the account #
# system.                                                              #
########################################################################

if (must.set) BalObjV <- set.BalObj(BalObjV, BalConstr)

# Default: cumulate variances on *all* constraints
n.constr <- nrow(BalConstr$obj.pos.mat)
if (is.null(i.constr)) {
     i.constr <- 1:n.constr
    }
# Map operator codes into appropriate R functions
op.fun <- function(op.code, obj) switch(op.code, `1`=cbind(rowSums(obj)), `2`=rbind(colSums(obj)), `3`=matrix(sum(obj)), `4`=obj)
# NOTE: it seems one could substitute *all* the matricial structures of BalConstr
#       (which involve 0s & NAs) with lists: it would be easier and lighter...
cumvar <- lapply(i.constr, function(i) {
                             obj.pos <- (BalConstr$obj.pos.mat[i, ] > 0)
                             objs <- BalConstr$obj.pos.mat[i, obj.pos]
                             sgns <- BalConstr$sgn.mat[i, obj.pos]
                             op.codes <- BalConstr$op.code.mat[i, obj.pos]
                             n.objs <- length(objs)
                             res <- Reduce(`+`, lapply(1:n.objs, function(j.obj)
                                                                          op.fun(op.codes[j.obj], BalObjV[[objs[j.obj]]])))
                            })
cumvar
}


`Balance` <- function(BalConstr, BalObj, BalObjV, tol = 1E-7,
                      maxit = 5000, check.conv = 5, force = TRUE, verbose = TRUE, time = FALSE, ...){
##################################################################
# Solve the balancing system via Conjugate Gradient Algorithm.   #
# The weighted least squares solution is seeked.                 #
# NOTE: The balanced object list is returned.                    #
# NOTE: Why not providing the inverse variance weighted RSS too? #
#       Could turn out to be useful, e.g. for assessing the      #
#       statistical quality of the result.                       #
# NOTE: Argument 'time' and related elaborations are here only   #
#       temporarily and will be removed when Speakeasy Task      #
#       Force is over.                                           #
##################################################################

# If requested start timing
if (time) Tstart <- proc.time()

# Must write sanity checks on input parameters: e.g. tol > 0, maxit >= 1, 1 <= check.conv < maxit, ...
if ( (maxit < 1) || (check.conv < 1) || (check.conv > maxit)){
     stop("Arguments maxit and check.conv must satisfy: [maxit >= 1] and [1 <= check.conv < maxit]\n\n")
    }

BalObj  <- set.BalObj(BalObj,  BalConstr)
BalObjV <- set.BalObj(BalObjV, BalConstr)

# Check that objects actually need to be modified
misfit <- ConstrMisfit(BalConstr, BalObj, must.set = FALSE)
ok.constr <- sapply(misfit, function(m) max(abs(m)) < tol)
if (all(ok.constr)){
     cat("All balancing constraints (", nrow(BalConstr$obj.pos.mat),") fulfilled (at tolerance level tol = ",
            tol, ").\n\n",sep="")
     return(invisible(NULL))
    }

# Check that one can actually modify what needs to be
cumvar <- ConstrCumVar(BalConstr, BalObjV, must.set = FALSE)
ko.var.constr <- mapply(function(m, v) any( (abs(m) > tol) & (v <= 0) ), misfit, cumvar)
if (any(ko.var.constr)){
     Text.ko.var.constr <- if (BalConstr$skip.decl) BalConstr$true.v.constr[ko.var.constr] else BalConstr$v.constr[ko.var.constr]
     warning("Estimates with zero *cumulative* variance found in the following failed (macro) constraints:\n  ",
          paste(Text.ko.var.constr, collapse = "\n  "), 
          "\n\n  REMARK: Balance might not converge and residual discrepancies will exist for balanced estimates!\n\n", 
          immediate. = TRUE)
    }

# If needed, summarize account balancing system
if (verbose){
     # NOTE: can skip *all* Bal.input stuff, as:
       # dim.obj.constr *is pointless*
     true.misfit <- ConstrMisfit(BalConstr, BalObj, BalConstr$true.constr.i, must.set = FALSE)
     n.micro.eq <- Reduce(`+`, lapply(true.misfit, function(e) prod(dim(e))))
       # obj.eq.dim *is pointless*
     n.atom.est <- Reduce(`+`, lapply(BalObj, function(e) prod(dim(e))))
     n.ko.constr <- length(misfit) - sum(ok.constr)

     atom.misfit <- abs(unlist(true.misfit))[abs(unlist(true.misfit)) > 0]
     n.atom.misfit <- length(atom.misfit)          

     cat("- Number of broken *macro*  constraints (at tolerance level tol = ", tol, "): ", n.ko.constr, "\n", sep="")
     cat("- Number of broken *atomic* constraints (at tolerance level tol = ", tol, "): ", n.atom.misfit, "\n", sep="")
     cat("- Summary of constraint violation errors (absolute values):\n\n")
         print(summary(atom.misfit))
     cat("\n")

     cat("- Number of atomic constraints in the balancing system: ", n.micro.eq, "\n", sep="")
     cat("- Number of atomic estimates to be balanced:            ", n.atom.est, "\n\n", sep="")

     vsep <- paste(rep("-", 0.12 * getOption("width")), collapse = "")
     cat("\n", vsep, " Conjugate Gradient Loop ", vsep, "\n\n", sep="")
    }

# Initialize Conjugate Gradient loop #
wpc <- lapply(cumvar, function(e) {
                             out <- 1/sqrt(e)
                             # as this can only occur for 0 misfits
                             out[!is.finite(out)] <- 0
                             out})
phi <- rho <- mapply(`*`, misfit, wpc, SIMPLIFY = FALSE)
lam <- lapply(misfit, function(e) e*0)

w.res <- unlist(phi)

n.objs <- nrow(BalConstr$constr.pos.obj)
i.obj <- 1:n.objs

# Start Conjugate Gradient loop #
iter <- check.iter <- 0
counter <- check.conv - 1  # check on a check.conv-fold basis + last
oldcontr <- Inf
osc <- FALSE
repeat ({
       iter <- iter + 1
       counter <- counter + 1
       tout <- lapply(i.obj, function(i) {
                      eqs <- BalConstr$constr.pos.obj[i, ]
                      j.eq <- which(eqs > 0)
                      tout.i <- lapply(j.eq, function (j, i) {
                                       s <- BalConstr$sgn.pos.obj[i, j]
                                       p <- phi[[eqs[j]]]
                                       w <- wpc[[eqs[j]]]
                                       return(s * p * w)
                                    }, i = i)
                      Reduce(`hwSum`, tout.i)
                    })

       # Compute misfit for weighted tout  (NOTE: reuse symbols)
       tout <- mapply(`*`, tout, BalObjV, SIMPLIFY = FALSE)
       misfit <- ConstrMisfit(BalConstr, tout, must.set = FALSE)

       A.phi <- mapply(`*`, misfit, wpc, SIMPLIFY = FALSE)
       rho.ss <- Reduce(`+`, lapply(rho, function(e) sum(e^2)))
       alpha <- rho.ss / Reduce(`sum`, mapply(`*`, phi, A.phi, SIMPLIFY = FALSE))

       rho <- mapply(function(a, b, c) a - b*c, rho, alpha, A.phi, SIMPLIFY = FALSE)
       lam <- mapply(function(a, b, c) a + b*c, lam, alpha, phi, SIMPLIFY = FALSE)

       rho.ss.next <- Reduce(`+`, lapply(rho, function(e) sum(e^2)))

       beta <- rho.ss.next/rho.ss
       phi <- mapply(function(a, b, c) a + b*c, rho, beta, phi, SIMPLIFY = FALSE)

       if ((counter == check.conv) || (iter >= maxit)) {
           check.iter <- check.iter + 1
           tout <- lapply(i.obj, function(i) {
                      eqs <- BalConstr$constr.pos.obj[i, ]
                      j.eq <- which(eqs > 0)
                      tout.i <- lapply(j.eq, function (j, i) {
                                       s <- BalConstr$sgn.pos.obj[i, j]
                                       l <- lam[[eqs[j]]]
                                       w <- wpc[[eqs[j]]]
                                       return(s * l * w)
                                    }, i = i)
                      Reduce(`hwSum`, tout.i)
                    })

           # Compute misfit for weighted tout  (NOTE: reuse symbols)
           tout <- mapply(`*`, tout, BalObjV, SIMPLIFY = FALSE)
           misfit <- ConstrMisfit(BalConstr, tout, must.set = FALSE)
           w.res.cur <- unlist(mapply(`*`, misfit, wpc, SIMPLIFY = FALSE))

           # Re-set check counter
           counter <- 0

           contr <- sum(abs(w.res - w.res.cur))   # ORIGINAL stopping criterion ## NOTE (12/04/2018): Perhaps better for DBE applications ##
           # NEW stopping criterion #
               # I think using max(abs()) is by far better:
               # 1) makes 'tol' more understandable to users setting its value
               # 2) makes Conjugate Gradient less prone to numerical instability,
               #    especially when the number of atomic estimates to be balanced
               #    gets huge.
           # contr <- max(abs(w.res - w.res.cur))
           msg <- paste("Check: ", check.iter, " - iteration: ", iter, " [residual error: ", contr, "]", sep="")

           if (contr > oldcontr) {
               msg <- paste(msg, " *system is oscillating!*", sep="")
               osc <- TRUE
            }

           if (verbose) cat(paste(msg, "\n", sep=""))

           if (contr < tol) {
               if (verbose) cat(paste("\n*Convergence achieved* in ", iter, " iterations [residual error: ", contr, "]\n\n", sep=""))
               break
            }

           if (iter >= maxit) {
               ko.msg <- paste("\n*Failed to converge* in ", iter, " iterations [residual error: ", contr, "]\n\n", sep="")
               if (force && !osc) {
                   warning(ko.msg)
                   break
                }
               else stop(ko.msg)
            }
            
           oldcontr <- contr
           osc <- FALSE
        }
    })
# Conjugate Gradient loop End #

BalObjQ <- BalObj
BalObjQ[] <- mapply(`-`, BalObj, tout, SIMPLIFY = FALSE)


# ##### WEIGHTED Residual Sum of Squares ##### START
# s2 <- lapply(BalObjV, function(e) {
                             # out <- 1/e
                             # # as this can only occur for 0 misfits
                             # out[!is.finite(out)] <- 0
                             # out})
# s <- lapply(BalObjV, function(e) {
                             # out <- 1/sqrt(e)
                             # # as this can only occur for 0 misfits
                             # out[!is.finite(out)] <- 0
                             # out})
# cat("RSS    = ", Reduce(`+`, mapply(function(x, y) sum((x-y)^2), BalObj, BalObjQ, SIMPLIFY = FALSE)), "\n")
# cat("RSS.s  = ", Reduce(`+`, mapply(function(x, y, w) sum(w*(x-y)^2), BalObj, BalObjQ, s, SIMPLIFY = FALSE)), "\n")
# cat("RSS.s2 = ", Reduce(`+`, mapply(function(x, y, w) sum(w*(x-y)^2), BalObj, BalObjQ, s2, SIMPLIFY = FALSE)), "\n\n")
# ##### WEIGHTED Residual Sum of Squares ##### END


attr(BalObjQ, "o.names")[] <- paste("Q", attr(BalObjQ, "o.names"), sep="")
BalObjQ <- reset.BalObj(BalObjQ)

# If requested stop timing
if (time) {
     Tstop <- proc.time()
     attr(BalObjQ, "time") <- (Tstop - Tstart)
    }

BalObjQ
}


`write.BalObj` <- function(BalObj, path = if (interactive()) choose.dir(getwd(), caption = "SELECT OUTPUT FOLDER"),
                           file.ext = "csv", sep = ";", dec = '.',
                           verbose = FALSE, ...){
#######################################################################
# Export balanced (or to-be-balanced) objects (i.e. numeric matrices) #
# stored inside BalObj list. Each object is written to a file whose   #
# name exactly matches the object name.                               #
#######################################################################

  obj.names <- names(BalObj)
  for (i.obj in seq_along(BalObj)){ 
       write.table(x = BalObj[[i.obj]],
                   file = file.path(path, paste(obj.names[i.obj], file.ext, sep=".")),
                   quote = FALSE,
                   sep = sep,
                   dec = dec,
                   row.names = FALSE,
                   col.names = FALSE,
                   ...)
    }

  if (verbose) {
      cat("- Number of objects written: ", length(BalObj), "\n\n")
    }

  invisible(NULL)
}


`hwSum` <- function(a, b){
##########################################################
# Sum two matrix objects according to Speakeasy language #
# "high-wide convention".                                #
##########################################################
  # Both a and b must be matrices
  if (is.matrix(a) && is.matrix(b)){
      if (identical(dim(a), dim(b))) return(a + b)
      a1 <- (1 %in% dim(a))
      b1 <- (1 %in% dim(b))
      # 1) Either a or b are (row/col) vectors, but NOT both
      if (xor(a1, b1)){
          m <- if (a1) b else a
          v <- if (a1) a else b
          is.row.v <- (ncol(v) > 1)
          is.scalar.v <- (prod(dim(v)) == 1)      # Copes with a(m,n) + b(1,1), where m > 1 and n > 1 (Francesco's code is flawed here)
          is.conform.v <- if (is.row.v) {ncol(v)==ncol(m)} else {nrow(v)==nrow(m)}
          if (is.conform.v || is.scalar.v){
              out <- m + matrix(v, nrow=nrow(m), ncol=ncol(m), byrow=is.row.v)
              return(out)
            }
          else {
              return(a + b)          
            }
        }
      # 2) Both a and b are (row/col) vectors
      if (a1 && b1){
          if (any(xor( (dim(a)==1), (dim(b)==1) ))){
              if (prod(dim(a))==1) return(b + a[,])
              if (prod(dim(b))==1) return(a + b[,])
              v.col <- if (nrow(a) > 1) a else b
              v.row <- if (ncol(a) > 1) a else b
              out <- outer(as.numeric(v.col), as.numeric(v.row), "+")
              return(out)
            }
          else {
              return(a + b)
            }
        }
      # 3) Neither a nor b are (row/col) vectors
      return(a + b)
    }
  else {
      stop("Both input objects (a and b) must be matrices!\n\n")
    }
}


`check.ZeroVarConstr` <- function(BalConstr, BalObj, BalObjV, tol = 1E-7, verbose = TRUE, ...){
###########################################################################
# In case broken constraints exist that involve only unalterable objects, #
# i.e. zero variance estimates, report on them. These constraints cannot  #
# possibly be fixed by the balancing algorithm!                           #
#                                                                         #
# NOTE: The array indices of all atomic constraints involved in detected  #
#       pathological macro constraints are returned.                      #
#                                                                         #
# NOTE: This pathology is signaled by a warning massage when executing    #
#       main function Balance().                                          #
###########################################################################

BalObj  <- set.BalObj(BalObj,  BalConstr)
BalObjV <- set.BalObj(BalObjV, BalConstr)

# Check that objects actually need to be modified
misfit <- ConstrMisfit(BalConstr, BalObj, must.set = FALSE)
ok.constr <- sapply(misfit, function(m) max(abs(m)) < tol)
if (all(ok.constr)){
     cat("All balancing constraints (", nrow(BalConstr$obj.pos.mat),") fulfilled (at tolerance level tol = ",
            tol, ").\n\n",sep="")
     return(invisible(NULL))
    }

# Check if one can actually modify what needs to be
cumvar <- ConstrCumVar(BalConstr, BalObjV, must.set = FALSE)

ko.var.constr <- mapply(function(m, v) any( (abs(m) > tol) & (v <= 0) ), misfit, cumvar)

if (!any(ko.var.constr)){
     cat("No failed balancing constraints involve zero *cumulative* variance estimates.\n\n",sep="")
     return(invisible(NULL))
    } else {
     Text.ko.var.constr <- if (BalConstr$skip.decl) BalConstr$true.v.constr[ko.var.constr] else BalConstr$v.constr[ko.var.constr]
     cat("Estimates with zero *cumulative* variance found in the following failed (macro) constraints:\n  ",
         paste(Text.ko.var.constr, collapse = "\n  "), 
         "\n\nREMARK: Residual discrepancies will exist for balanced estimates or balancing algorithm might not converge!\n\n", sep = "")

     atom.ko.var.constr <- mapply(function(m, v) which( (abs(m) > tol) & (v <= 0) , arr.ind = TRUE),
                                  misfit[ko.var.constr], cumvar[ko.var.constr], SIMPLIFY = FALSE)

     if (verbose) {
         cat("Array indices of involed atomic constraints:\n\n")
         print(atom.ko.var.constr)
        }

     return(invisible(atom.ko.var.constr))
    }
}


`BalanceR` <- function(constr.file = if (interactive()) choose.files(multi = FALSE, caption = "SELECT CONSTRAINTS FILE"),
                       obj.dir     = if (interactive()) choose.dir(getwd(), caption = "SELECT OBJECTS FOLDER"),
                       obj.var.dir = if (interactive()) choose.dir(getwd(), caption = "SELECT VARIANCES FOLDER"),
                       out.obj.dir = if (interactive()) choose.dir(getwd(), caption = "SELECT OUTPUT FOLDER"),
                       tol = 1E-7, maxit = 5000, check.conv = 5, force = TRUE,
                       file.ext = "csv", sep = ";", dec = '.',
                       skip.decl = TRUE, verbose = TRUE, time = FALSE, ...){
##########################################################################
# A driver running a balancing procedure as a whole.                     #
# NOTE: When out.obj.dir is NULL (the default), a list storing balanced  #
#       output objects is returned. If, instead, the path to an existing #
#       directory is provided, balanced output objects are exported to   #
#       that directory.                                                  #
# NOTE: Argument 'time' and related elaborations are here only           #
#       temporarily and will be removed when Speakeasy Task Force is     #
#       over.                                                            #
##########################################################################

BalConstr <- read.BalConstr(constr.file, skip.decl = skip.decl, verbose = verbose, time = time)
BalObj <- read.BalObj(obj.dir, file.ext = file.ext, sep = sep, dec = dec)
BalObjV <- read.BalObj(obj.var.dir, file.ext = file.ext, sep = sep, dec = dec)

dont.write <- ( is.null(out.obj.dir) || is.na(out.obj.dir) )

BalObjQ <- Balance(BalConstr, BalObj, BalObjV,
                   tol = tol, maxit = maxit, check.conv = check.conv, force = force, verbose = verbose, time = time)

if (verbose & time & !is.null(Tconstr <- attr(BalConstr, "time"))) {
     cat("\n*Elaboration time*\n")
     print(ttime <- Tconstr + attr(BalObjQ, "time"))
     cat("\n\n")
}

attr(BalObjQ, "time") <- NULL

if (dont.write) {
     return(BalObjQ)
    }
else {
     write.BalObj(BalObj = BalObjQ, path = out.obj.dir, file.ext = file.ext, sep = sep, dec = dec)
    }
}
