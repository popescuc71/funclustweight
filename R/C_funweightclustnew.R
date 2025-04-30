#********************************* REQUIRED LIBRARIES *********************************#

# library(stats)
# #library(funHDDC)
# library(fda)
# library(mclust)
# library(tclust)
# library(stringr)
# library(MASS)

#********************************* Clustering Function *********************************#
# /* BEGIN FROM FUNHDDC (MODIFIED) */
# /*
# * Authors: Bouveyron, C. Jacques, J.
# * Date Taken: 2022-01-01
# * Original Source: funHDDC (modified)
# * Address: https://github.com/cran/funHDDC
# *
# */
funclustweight  <-
  function(datax,datay, K=1:10, model="AkjBkQkDk", modely="VVV",known=NULL, threshold=0.1, itermax=200, eps=1e-6, init='random',
           criterion="bic", d_select="Cattell", init.vector=NULL,
           show=TRUE, mini.nb=c(5, 10), min.individuals=2, mc.cores=1, nb.rep=2,
           keepAllRes=TRUE, kmeans.control = list(), d_max=100, d_range=2,cmtol=1e-10,
           cmmax=10,verbose = TRUE){
    

    #Options removed from call
    com_dim <- NULL

    noise.ctrl <- 1e-8

    #
    # CONTROLS
    #

    call <- match.call()  # macthes values to prameters

    .reg_hddc_control(call)
    # Control of match.args:
    criterion <- .reg_myAlerts( criterion, "criterion", "singleCharacterMatch.arg",
                            "funclustweight: ", c("bic", "icl") )

    d_select <- .reg_myAlerts( d_select, "d_select", "singleCharacterMatch.arg",
                           "funclustweight: ", c("cattell", "bic", "grid") )
    init <- .reg_myAlerts( init, "init", "singleCharacterMatch.arg",
                       "funclustweight: ", c('random', 'kmeans',
                                   'mini-em', "vector") )
    # We get the model names, properly ordered
    model <- .reg_hdc_getTheModel(model, all2models = TRUE)
    modely<-.reg_hdc_getTheModely(modely, all2models = TRUE)
#@@@@@@@@@@@@need models for y
    #nb.rep parameter
    if ( (init == "random") & (nb.rep < 20) ){
      nb.rep <- 20
    }

    # set verbose to false if using multiple cores, no point in timing
    if(mc.cores > 1) {
      verbose <- FALSE
    }

    # kmeans controls
   
    kmeans.control <- .reg_default_kmeans_control(kmeans.control)


    BIC <- ICL <- c()
    fdobj = datax  #for univariate  model it is NOT a list
    fdobjy = datay
    Wlist<-list()
    if (!inherits(fdobj, 'list'))
    {
      x <- t(fdobj$coefs)
      p <- ncol(x)

      #Constructing the matrix with inner products of the basis functions
      W <- inprod(fdobj$basis,fdobj$basis)
      W[W < 1e-15] <- 0
      #Construction of the triangular matrix of Choleski
      W_m <- chol(W)
      dety=det(W)
      Wlist<-list(W=W,
                  W_m=W_m,
                  dety=dety
      )
    }
    else {
      x = t(fdobj[[1]]$coefs);
      for (i in 2:length(fdobj)) x <- cbind(x,t(fdobj[[i]]$coefs))
      p <- ncol(x)
      for ( i in 1:length(fdobj) ){
        name <- paste('W_var', i, sep = '')
        #Constructing the matrix with inner products of the basis functions
        W_fdobj <- inprod(fdobj[[i]]$basis, fdobj[[i]]$basis)
        assign(name, W_fdobj)
      }

      #Add 0 ? left and right of  W before constructing the matrix phi
      prow <- dim(W_fdobj)[[1]]
      pcol <- length(fdobj) * prow
      W1 <- cbind( W_fdobj,
                   matrix( 0, nrow = prow, ncol = ( pcol - ncol(W_fdobj) ) ) )
      W_list <- list()
      for ( i in 2:( length(fdobj) ) ){
        W2 <- cbind( matrix( 0, nrow = prow, ncol = (i - 1) * ncol(W_fdobj) ),
                     get( paste('W_var', i, sep = '') ),
                     matrix( 0, nrow = prow, ncol = ( pcol - i * ncol(W_fdobj) ) ) )
        W_list[[i - 1]] <- W2
      }

      #Constructing the matrix phi
      W_tot <- rbind(W1,W_list[[1]])
      if (length(fdobj) > 2){
        for( i in 2:(length(fdobj) - 1) ){
          W_tot <- rbind(W_tot, W_list[[i]])
        }
      }
      W_tot[W_tot < 1e-15] <- 0
      #Construction of the triangular matrix of Choleski
      W_m <- chol(W_tot)
      dety<-det(W_tot)
      Wlist<-list(W=W_tot,
                  W_m=W_m,
                  dety=dety
      )
      }


    #
    # Preparing the parallel
    #

    # bic and grid dimension selection do not bother with a threshold
    if(d_select == "bic"){
      threshold <- "bic"
    } else if(d_select == "grid"){
      threshold <- "grid"
    }

    if(max( table(K) ) > 1) warning("The number of clusters, K, is made unique (repeated values are not tolerated).")
    K <- sort( unique(K) )

    if(d_select == "grid") {
      # set up the grid of values for d_set
      if(length(K) > 1) stop("funclustweight: Using d_select='grid' K must be only 1
                             value (not a list).")
      mkt_list <- list(model = model, modely=modely, K = K, threshold = threshold)
      for(i in 1:K) mkt_list[[paste0("d", i)]] <- d_range
      mkt_Expand <- do.call(expand.grid, mkt_list)

    } else {
      mkt_list <- list(model = model, modely = modely, K = K, threshold = threshold)
      # default needs to be a number set as a number just in case
      for( i in 1:max(K) ) mkt_list[[paste0("d", i)]] <- 2
      mkt_Expand <- do.call(expand.grid, mkt_list)
    }
    mkt_Expand <- do.call( rbind, replicate(nb.rep, mkt_Expand,
                                            simplify = FALSE) )
    mkt_Expand <- rbind(c(), mkt_Expand) #no need for several runs for K==1

    model <- as.character(mkt_Expand$model)
    modely <- as.character(mkt_Expand$modely)
    K <- mkt_Expand$K
    threshold <- mkt_Expand$threshold
    d <- list()
    for( i in 1:max(K) ) {
      d[[i]] <- mkt_Expand[, 3 + i]
    }

    # timing should be initialized and total time details should be created
    if(verbose) {

      backspace.amount <- .reg_estimateTime("init")[2]

      mkt_Expand <- cbind( mkt_Expand, modelcount = 0:(nrow(mkt_Expand)-1) )
    }
    # We transform it into an univariate form
    mkt_univariate <- apply(mkt_Expand, 1, paste, collapse = "_")


    # Mon 'caller' for LTBM that will be used in multi-cores lapply
    hddcWrapper <- function(mkt_univariate, verbose,
                            start.time=0, totmod=1, backspace.amount=0, ...){


      mkt_splitted <- strsplit(mkt_univariate, "_")

      # we find model, K and threshold
      model <- sapply(mkt_splitted, function(x) x[1])
      modely <- sapply(mkt_splitted, function(x) x[2])
      K <- sapply( mkt_splitted, function(x) as.numeric(x[3]) )
      threshold <- sapply(mkt_splitted, function(x) ifelse( is.numeric(x[4]),
                                                            as.numeric(x[4]),
                                                            x[4]) )
     

      # note d_set is ignored unless run in d_select = "grid"
      d_set <- rep(2, K)
      for(i in 1:K) d_set[i] <- sapply(mkt_splitted,
                                       function(x) as.numeric(x[4+i]))

      # (::: is needed for windows multicore)
      res <- "unknown error"



      # removed try wrapper
      try(res <- .Reg_funhddc_main1(model = model,modely=modely, K = K,
                                 threshold = threshold,
                                 d_set = d_set,...), silent = T)
      if(verbose) {
        # update the timing
          modelcount <- sapply( mkt_splitted, function(x) as.numeric(x[length(mkt_splitted[[1]])]) )
          rv <- .reg_estimateTime(modelcount, start.time, totmod, backspace = backspace.amount)
      }
      res
    }

    # We reset the number of cores to use
    nRuns <- length(mkt_univariate)
    if(nRuns < mc.cores) mc.cores <- nRuns

    # We swicth to the right number of cores + a warning if necessary
    max_nb_of_cores <- parallel::detectCores()
    if(mc.cores > max_nb_of_cores){
      warning("The argument mc.cores is greater than its maximun.
              \nmc.cores was set to ", max_nb_of_cores)
      mc.cores <- max_nb_of_cores
    }


    #
    # Parallel estimations
    #

    start.time <- Sys.time() # get a starting time
    if(mc.cores == 1){
      # If there is no need for parallel, we just use lapply // in order to have the same output

      # details for timing should be passed to lapply


      par.output <- lapply(mkt_univariate, hddcWrapper, fdobj = fdobj,fdobjy=fdobjy,Wlist=Wlist, known=known,
                           method = d_select, eps = eps,
                           init = init, init.vector = init.vector,
                           mini.nb = mini.nb, min.individuals = min.individuals,
                           noise.ctrl = noise.ctrl, com_dim = com_dim,
                           kmeans.control = kmeans.control, d_max = d_max,
                           start.time = start.time,cmtol=cmtol,
                           cmmax=cmmax, verbose = verbose,
                           backspace.amount = backspace.amount, totmod = nrow(mkt_Expand))
      #  here mkt_univariate is applied to to fdobj

    } else {

      # we use parLapply:

      ## create clusters
      cl <- parallel::makeCluster(mc.cores)

      parallel::clusterEvalQ( cl, library(funclustweight) )
      parallel::clusterExport( cl, ls( environment() ), envir = environment() )

      ## run the parallel
      par.output <- NULL
     
      try(par.output <- parallel::parLapplyLB( cl, mkt_univariate, hddcWrapper,
                                             fdobj=fdobj,fdobjy=fdobjy,
                                             Wlist=Wlist, method=d_select,
                                             known=known,
                                             itermax=itermax,
                                             eps=eps, init=init,
                                             init.vector=ifelse(missing(init.vector),
                                                                NA, init.vector),
                                             mini.nb=mini.nb,
                                             min.individuals=min.individuals,
                                             noise.ctrl=noise.ctrl,
                                             com_dim=com_dim,
                                             kmeans.control=kmeans.control,
                                             d_max=d_max, start.time = start.time,
                                             cmtol=cmtol,
                                             cmmax=cmmax, verbose = FALSE) )

      ## Stop the clusters
      parallel::stopCluster(cl)

      if( is.null(par.output) ) stop("Unknown error in the parallel computing.
                                     Try mc.cores=1 to detect the problem.")

    }
    # The results are retrieved
    #

    getElement <- function(x, what, valueIfNull = -Inf){
      # atention if x is the null model
      if(length(x) == 1) return(valueIfNull)
      if(!is.list(x) && !what %in% names(x)) return(NA)
      x[[what]][length(x[[what]])]
    }

    getComment <- function(x){
      # we get the error message
      if(length(x) == 1) return(x)
      return("")
    }

    # All likelihoods
    LL_all <- sapply(par.output, getElement, what = "loglik")
    comment_all <- sapply(par.output, getComment)

    # If no model is valid => problem
    if( all( !is.finite(LL_all) ) ){
      warning("All models diverged.")
      allCriteria <- data.frame(model = model,modely= modely, K = K, threshold = threshold,
                                LL = LL_all, BIC = NA, comment = comment_all)
      res <- list()
      res$allCriteria <- allCriteria
      return(res)
    }

    # We select, for each (Q,K), the best run
    n <- nrow(mkt_Expand)
    modelKeep <- sapply(unique(mkt_univariate),
                        function(x) (1:n)[mkt_univariate==x][which.max(LL_all[mkt_univariate==x])])
    # => we select only the best models
    LL_all <- LL_all[modelKeep]
    comment_all <- comment_all[modelKeep]
    par.output <- par.output[modelKeep]
    BIC <- sapply(par.output, getElement, what = "BIC")
    ICL <- sapply(par.output, getElement, what = "ICL")
    comp_all <- sapply(par.output, getElement, what = "complexity",
                       valueIfNull = NA)
    model <- model[modelKeep]
    modely <- modely[modelKeep]
    threshold <- threshold[modelKeep]
    K <- K[modelKeep]
    d_keep <- list()
    for(i in 1:max(K)) {
      d_keep[[i]] <- d[[i]][modelKeep]
    }

    # We define the criterion of model selection
    CRIT <- switch(criterion,
                   bic = BIC,
                   icl = ICL)

    # The order of the results
    myOrder <- order(CRIT, decreasing = TRUE)

    # we keep the good model + creation of output
    qui <- which.max(CRIT)

    prms <- par.output[[qui]]
    prms$criterion <- CRIT[qui]
    names(prms$criterion) <- criterion

    # Other output
    prms$call <- call
    # We add the complexity
    names(comp_all) <- mkt_univariate[modelKeep]
    prms$complexity_allModels <- comp_all

    # Display
    if(show){
      if(n > 1) cat("funclustweight: \n")

      model2print <- sapply(model,
                            function(x) sprintf( "%*s", max( nchar(model) ), x) )
      modely2print <- sapply(modely,
                             function(x) sprintf( "%*s", max( nchar(model) ), x) )
      K2print <- as.character(K)
      K2print <- sapply(K2print,
                        function(x) sprintf( "%*s", max( nchar(K2print) ), x) )
      thresh2print <- as.character(threshold)
      thresh_width <- max( nchar(thresh2print) )
      thresh2print <- sapply(thresh2print,
                             function(x) sprintf( "%s%s", x,
                                                  paste0( rep("0",
                                                             thresh_width - nchar(x) ),
                                                         collapse = "") ) )

      # we make a data.frame
      d_names <- c()
      myResMat <- cbind( model2print[myOrder], modely2print[myOrder], K2print[myOrder],
                         thresh2print[myOrder],
                         .reg_addCommas(prms$complexity_allModels[myOrder]),
                         .reg_addCommas(CRIT[myOrder]) )
      if(d_select == "grid") {
        # only display d_set values when they are relevant
        for(i in 1:max(K)) {
          d2print <- as.character(d_keep[[i]])
          d2print <- sapply( d2print,
                             function(x) sprintf( "%*s", max( nchar(d2print) ), x) )
          myResMat <- cbind(myResMat, d2print[myOrder])
          d_names <- append( d_names, paste0("d", i) )
        }
      }
      myResMat <- cbind(myResMat, comment_all[myOrder])

      myResMat <- as.data.frame(myResMat)
      names(myResMat) <- c("model", "modely", "K", "threshold", "complexity",
                           toupper(criterion), d_names, "comment")
      row.names(myResMat) <- 1:nrow(myResMat)

      # if no problem => no comment
      if( all(comment_all == "") ) myResMat$comment <- NULL

      print(myResMat)

      msg <- switch(criterion, bic = "BIC", icl = "ICL")
      cat("\nSELECTED: model ", prms$model, "-", prms$modely, " with ", prms$K, " clusters.\n")
      cat("Selection Criterion: ", msg, ".\n", sep="")

    } # end of if(show)

    # We also add the matrix of all criteria
    allCriteria <- data.frame( model = model[myOrder],
                               modely = modely[myOrder],
                               K = K[myOrder],
                               threshold = threshold[myOrder],
                               LL = LL_all[myOrder],
                               complexity = prms$complexity_allModels[myOrder],
                               BIC = BIC[myOrder],
                               ICL = ICL[myOrder],
                               rank = 1:length(myOrder) )

    # we add the comments if necessary
    if( any(comment_all != "") ) allCriteria$comment <- comment_all[myOrder]
    prms$allCriteria <- allCriteria

    # If all the results are kept
    if(keepAllRes){
      all_results <- par.output
      names(all_results) <- mkt_univariate[modelKeep]
      prms$all_results <- all_results
    }

    # Other stuff

    prms$threshold <- threshold[qui]
    prms$datax<-fdobj
    prms$datay<-fdobjy
    class(prms)<-"t-funHDDC"
    return(prms)
  } # end of tfunHDDC

.Reg_funhddc_main1 <- function(fdobj,fdobjy,Wlist, K,model,modely, itermax=200,threshold,
                            method, eps=1e-6, init, init.vector, mini.nb,
                            min.individuals, noise.ctrl, com_dim=NULL,
                            kmeans.control, d_max, d_set=c(2,2,2,2),known=NULL,cmtol=cmtol,
                            cmmax=cmmax, ...){

  ModelNames <- c("AKJBKQKDK","AKJBQKDK", "AKBKQKDK", "ABKQKDK", "AKBQKDK", "ABQKDK")
  ModelNamesy <- c("VVV")

  
  if (!inherits(fdobj, 'list')) {DATA <- t(fdobj$coefs)}
  else {DATA = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) DATA <- cbind(DATA,t(fdobj[[i]]$coefs))}
  if (!inherits(fdobjy, 'list')) { DATAy <- t(fdobjy$coefs) # THIS IS the univariate CASE
  } else {DATAy = t(fdobjy[[1]]$coefs); for (i in 2:length(fdobjy))
    DATAy = cbind(DATAy,t(fdobjy[[i]]$coefs))}


  p <- ncol(DATA)
  N <- nrow(DATA)
  q<-ncol(DATAy)
  bigDATA <-rbind(Wlist$W%*%t(DATA), rep(1,N))#it is (p+1)XN
  com_ev <- NULL

  # We set d_max to a proper value
  d_max <- min(N, p, d_max)
  #################################clasification######################
  ########################################################
  # /*
  # * Authors: Andrews, J. Wickins, J. Boers, N. McNicholas, P.
  # * Date Taken: 2023-01-01
  # * Original Source: teigen (modified)
  # * Comment: Known and training are treated as in teigen, matchtab and matchit
  # *          are also treated like in teigen
  # * Address: https://github.com/cran/teigen
  # *
  # */
  if(is.null(known)){
    #this is the  clustering case
    clas <- 0
    kno <- NULL
    testindex <- NULL
  }
  else
    #various sorts of classification
    {
    if(length(known)!=N){
      return("Known classifications vector not given, or not the same length as the number of samples (see help file)")
    }else
    {

      if(!anyNA(known)){
        #Gs <- length(unique(known))
        warning("No NAs in 'known' vector supplied, all values have known classification (parameter estimation only)")
        testindex <- 1:N
        kno <- rep(1, N)
        unkno <- (kno-1)*(-1)
        K <- length(unique(known))
        init.vector <- as.numeric(known)
        init <-"vector"
      }
      else{
        training <- which(!is.na(known))
        testindex <- training
        kno <- vector(mode="numeric", length=N)
        kno[testindex] <- 1
        unkno <- (kno-1)*(-1)##it is 0 for training and 1 for the rest
      }
    }

    clas <- 1

  }
  ############################################
  #################################################

  if (K > 1){
    t <- matrix(0, N, K)
    tw <- matrix(0, N, K)
    ##############################################
    #initialization for t
    if(init == "vector"){
      if(clas>0){
        cn1=length(unique(known[testindex]))
        matchtab=matrix(-1,K,K)
        matchtab[1:cn1,1:K] <- table(known,init.vector)
        rownames(matchtab)<-c(rownames(table(known,init.vector)),which(!(1:K %in% unique( known[testindex]))))
        matchit<-rep(0,1,K)
        while(max(matchtab)>0){

          ij<-as.integer(which(matchtab==max(matchtab), arr.ind=T)[1,2])
          ik<-which.max(matchtab[,ij])
          matchit[ij]<-as.integer(rownames(as.matrix(ik)))
          matchtab[,ij]<-rep(-1,1,K)
          matchtab[as.integer(ik),]<-rep(-1,1,K)

        }
        matchit[which(matchit==0)]=which(!(1:K %in% unique(matchit)))
          initnew <- init.vector
          for(i in 1:K){
            initnew[init.vector==i] <- matchit[i]
          }
          init.vector <- initnew
        }


      for (i in 1:K) t[which(init.vector == i), i] <- 1
    }
    else {if (init == "kmeans") {
      kmc <- kmeans.control
      DATAbig <-cbind(DATA,DATAy)
      cluster <- kmeans(DATAbig, K, iter.max = kmc$iter.max, nstart = kmc$nstart,
                        algorithm = kmc$algorithm, trace = kmc$trace)$cluster
      if(clas>0)
        {
        cn1=length(unique(known[testindex]))
        matchtab=matrix(-1,K,K)
        matchtab[1:cn1,1:K] <- table(known,cluster)
        rownames(matchtab)<-c(rownames(table(known,cluster)),which(!(1:K %in% unique( known[testindex]))))
        matchit<-rep(0,1,K)
        while(max(matchtab)>0){

          ij<-as.integer(which(matchtab==max(matchtab), arr.ind=T)[1,2])
          ik<-which.max(matchtab[,ij])
          matchit[ij]<-as.integer(rownames(as.matrix(ik)))
          matchtab[,ij]<-rep(-1,1,K)
          matchtab[as.integer(ik),]<-rep(-1,1,K)

        }
        matchit[which(matchit==0)]=which(!(1:K %in% unique(matchit)))
          knew <- cluster
          for(i in 1:K){
            knew[cluster==i] <- matchit[i]
          }
          cluster <- knew

      }


      for (i in 1:K)
        t[which(cluster == i), i] <- 1

      #A@5 initialize t with the result of k-means
    } else
      {



     if (init=="mini-em"){

      prms_best <- 1
      for (i in 1:mini.nb[1]){
        prms <- .Reg_funhddc_main1(fdobj,fdobjy,Wlist, K, known=known, model = model,modely = modely,
                                threshold = threshold, method =  method,
                                itermax = mini.nb[2],
                                init.vector = 0,
                                init = 'random', mini.nb = mini.nb,
                                min.individuals = min.individuals,
                                noise.ctrl = noise.ctrl,
                                kmeans.control = kmeans.control,
                                com_dim = com_dim, d_max = d_max, d_set = d_set,cmtol=cmtol,
                                cmmax=cmmax)
        if(length(prms) != 1){
          if (length(prms_best) == 1) prms_best <- prms
          else if (prms_best$loglik[length(prms_best$loglik)] <
                   prms$loglik[length(prms$loglik)]) prms_best <- prms
        }
      }

      if (length(prms_best) == 1) return(1)
      t <- prms_best$posterior
      if(clas>0)
      {
        cluster <- max.col(t)
        cn1=length(unique(known[testindex]))
        matchtab=matrix(-1,K,K)
        matchtab[1:cn1,1:K] <- table(known,cluster)
        rownames(matchtab)<-c(rownames(table(known,cluster)),which(!(1:K %in% unique( known[testindex]))))
        matchit<-rep(0,1,K)
        while(max(matchtab)>0){

          ij<-as.integer(which(matchtab==max(matchtab), arr.ind=T)[1,2])
          ik<-which.max(matchtab[,ij])
          matchit[ij]<-as.integer(rownames(as.matrix(ik)))
          matchtab[,ij]<-rep(-1,1,K)
          matchtab[as.integer(ik),]<-rep(-1,1,K)

        }
        matchit[which(matchit==0)]=which(!(1:K %in% unique(matchit)))
          knew <- cluster
          for(i in 1:K){
            knew[cluster==i] <- matchit[i]
          }
          cluster <- knew
          for (i in 1:K)
            t[which(cluster == i), i] <- 1
        }

      }


     else {if (init=="random"){ #INIT IS RANDOM

      t <- t( rmultinom( N, 1, rep(1 / K, K) ) ) # some multinomial
      compteur <- 1
      while(min( colSums(t) ) < 1 && (compteur <- compteur + 1) < 5)
        t <- t( rmultinom( N, 1, rep(1 / K, K) ) )
      if(min( colSums(t) ) < 1)
        return("Random initialization failed (n too small)")
      if(clas>0)
      {
        cluster <- max.col(t)
        cn1=length(unique(known[testindex]))
        matchtab=matrix(-1,K,K)
        matchtab[1:cn1,1:K] <- table(known,cluster)
        rownames(matchtab)<-c(rownames(table(known,cluster)),which(!(1:K %in% unique( known[testindex]))))
        matchit<-rep(0,1,K)
        while(max(matchtab)>0){

          ij<-as.integer(which(matchtab==max(matchtab), arr.ind=T)[1,2])
          ik<-which.max(matchtab[,ij])
          matchit[ij]<-as.integer(rownames(as.matrix(ik)))
          matchtab[,ij]<-rep(-1,1,K)
          matchtab[as.integer(ik),]<-rep(-1,1,K)

        }
        matchit[which(matchit==0)]=which(!(1:K %in% unique(matchit)))
          knew <- cluster
          for(i in 1:K){
            knew[cluster==i] <- matchit[i]
          }
          cluster <- knew
          for (i in 1:K)
            t[which(cluster == i), i] <- 1
        }

      }

    }}}}
   else {
    t <- matrix(1, N, 1) # K IS 1

   }
  if(clas>0){
    t <- unkno*t
    for(i in 1:N){
      if(kno[i]==1){
        t[i, known[i]] <- 1
      }
    }
  }



  likely <- c()
  I <- 0
  test <- Inf
  while ( (I <- I + 1) <= itermax && test >= eps ){
    # loops here until itermax or test <eps
    # Error catching
    if (K > 1){
      if( any( is.na(t) ) ) return("unknown error: NA in t_ik")

      if( any(colSums(t > 1 / K) < min.individuals) )
        return("pop<min.individuals")
    }

    m <- .funclustweight_m_step(fdobj,bigDATA,fdobjy,Wlist, K, t, model,modely,
                          threshold, method, noise.ctrl, com_dim, d_max, d_set,cmtol=cmtol,
                          cmmax=cmmax) # mstep is applied here
    
    if(any(is.nan(m$icovy))) return("Model Error: icovy contains nans")


    to <- .funclustweight_e_step(fdobj,bigDATA,fdobjy,Wlist, m,clas, known, kno) # E-step is applied here
    L <- to$L #  this s the likelihood
    t <- to$t# this is the new t



    #likelihood contrib=L
    likely[I] <- L

    if (I == 2) test <- abs(likely[I] - likely[I - 1])
    else if (I > 2)
    {
      lal <- (likely[I] - likely[I - 1]) / (likely[I - 1] - likely[I - 2])
      lbl <- likely[I - 1] + (likely[I] - likely[I - 1]) / (1.0 - lal)
      test <- abs(lbl - likely[I - 1])
    }

  }



  # a
  if ( model%in%c("AKBKQKDK", "AKBQKDK", "AKBKQKD", "AKBQKD") ) {
    a <- matrix(m$a[, 1], 1, m$K, dimnames=list(c("Ak:"), 1:m$K))
  } else if(model=="AJBQD") {
    a <- matrix(m$a[1, ], 1, m$d[1], dimnames=list(c('Aj:'), paste('a', 1:m$d[1], sep='')))
  } else if ( model%in%c("ABKQKDK", "ABQKDK", "ABKQKD", "ABQKD", "ABQD") ) {
    a <- matrix(m$a[1], dimnames=list(c('A:'), c('')))
  } else a <- matrix(m$a, m$K, max(m$d), dimnames=list('Class'=1:m$K, paste('a', 1:max(m$d), sep='')))

  # b
  if ( model%in%c("AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBQKD", "AKBQKD", "ABQKD",
  "AJBQD", "ABQD") ) {
    b <- matrix(m$b[1], dimnames=list(c('B:'), c('')))
  } else b <- matrix(m$b, 1, m$K, dimnames=list(c("Bk:"), 1:m$K))

  # d, mu, prop
  d <- matrix(m$d, 1, m$K, dimnames=list(c('dim:'), "Intrinsic dimensions of the classes:"=1:m$K))
  mu <- matrix(m$mu, m$K, p, dimnames=list('Class'=1:m$K, 'Posterior group means:'=paste('V', 1:p, sep='')))
  prop <- matrix(m$prop, 1, m$K, dimnames=list(c(''), 'Posterior probabilities of groups'=1:m$K))

  # Other elements
  complexity <- .Reg.hdc_getComplexity(m, p,q)
  class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- 'hd' #A@5 alpha and eta will nedd class too
  cls <- max.col(t)#@ here the cluster is found

  ##
  converged = test < eps


  params <- list()
  params = c(params,list(
                Wlist=Wlist,
                model = model,
                modely = modely,
                K = K,
                d = d,
                a = a,
                b = b,
                mu = mu,
                prop = prop,
                ev = m$ev,
                Q = m$Q,
                Q1 = m$Q1,
                gam =m$gam,
                covy=m$covy,
                icovy=m$icovy,
                ldetcov=m$ldetcov,
                fpca = m$fpcaobj,
                loglik = likely[length(likely)],
                loglik_all = likely,
                posterior = t,
                class = cls,
                com_ev = com_ev,
                N = N,
                complexity = complexity,
                threshold = threshold,
                d_select = method,
                converged = converged))



  if(clas>0){
    params[["index"]] <- testindex
  }
  # We compute the BIC / ICL
  bic_icl <- .Reg_hdclassift_bic(params, p,q)
  params$BIC <- bic_icl$bic
  params$ICL <- bic_icl$icl

  # We set the class

  class(params) <- 'tfunHDDC'

  return(params)
} # end of .reg_funhddc_main1

#########################################################

######################################################################
############################################################
.funclustweight_e_step <- function(fdobj,bigDATA,fdobjy,Wlist, par,clas=0,known=NULL, kno=NULL){ # this is the E step
  if (!inherits(fdobj, 'list')) {x <- t(fdobj$coefs)}  # THIS IS the univariate CASE
  else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))} #NOT THIS
  if (!inherits(fdobjy, 'list')) { y <- t(fdobjy$coefs) # THIS IS the univariate CASE
  } else {y = t(fdobjy[[1]]$coefs); for (i in 2:length(fdobjy))
    y = cbind(y,t(fdobjy[[i]]$coefs))}
  bigx=t(bigDATA)
   p <- ncol(x)
   q<- ncol(y)
  N <- nrow(x)
  K <- par$K
  Ny <- nrow(y)
  py <- ncol(y)
  a <- par$a
  b <- par$b
  mu <- par$mu
  d <- par$d
  prop <- par$prop
  Q <- par$Q
  Q1 <- par$Q1
  b[b<1e-6] <- 1e-6
  icovy=par$icovy
  ldetcov=par$ldetcov
  gam=par$gam
  dety<-Wlist$dety
  if(clas>0){
    unkno <- (kno-1)*(-1)
  }
  ##################################################
  #########################################################
  ###########################################
  t <- matrix(0, N, K)

  mah_pen <- matrix(0, N, K)
  mah_pen1 <- matrix(0, N, K)
  # comp_pen1 <- matrix(0, N, K) # IAIN remove later
  K_pen <- matrix(0, K, N)
  num <- matrix(0, N, K)
  ft <- matrix(0, N, K)

  s <- rep(0, K)

  for (i in 1:K) {
    s[i] <- sum( log(a[i, 1:d[i]]) )

    Qk <- Q1[[i]]
    aki <- sqrt( diag( c( 1 / a[i, 1:d[i]], rep(1 / b[i], p - d[i]) ) ) )
    muki <- mu[i, ]

    Wki <- Wlist$W_m
    dety<-Wlist$dety

    mah_pen[, i] <- .reg_imahalanobis(x, muki, Wki ,Qk, aki)

    pqp=p+1
    # IAIN add again later
    res <- .C(".C_rmahalanobis", as.integer(N),as.integer(pqp),as.integer(q),
              as.integer(K),as.integer(i),
              as.double(bigx), as.double(y), as.double(gam[,,i]),
              as.double(icovy[,,i]),delta = rep(0,N), PACKAGE="funclustweight")
    # mah_pen1[, i] <- C_mahalanobis(bigx,y,gam[[i]], icovy[[i]]) # IAIN remove later
    mah_pen1[, i] <-res$delta
    K_pen[i, ]<--2*log(prop[i])+(p+q)*log(2*pi)+s[i] -log(dety)+ (p - d[i]) * log(b[i])+mah_pen[, i]+
                    +mah_pen1[, i]+ldetcov[i]

  }
  A <- -1/2*t(K_pen)
  L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))

  t <- matrix(0,N,K)
  for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))



  if(clas>0){
    t <- unkno*t
    for(i in 1:N){
      if(kno[i]==1){
        t[i, known[i]] <- 1
      }
    }
  }
  # if(any(is.nan(t))){
  #   break
  # }
  list(t = t,
       L = L)
} # end of .reg_funhddt_e_step1

.funclustweight_m_step  <- function(fdobj,bigDATA,fdobjy, Wlist, K, t, model,modely, threshold, method, noise.ctrl, com_dim, d_max, d_set,cmtol=cmtol,
                                   cmmax=cmmax){ #A@5 this is the M step
  if (!inherits(fdobj, 'list')) { x <- t(fdobj$coefs) # THIS IS the univariate CASE
  } else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
  # x is the coefficient in the fdobject
  if (!inherits(fdobjy, 'list')) { y <- t(fdobjy$coefs) # THIS IS the univariate CASE
  } else {y = t(fdobjy[[1]]$coefs); for (i in 2:length(fdobjy))
    y = cbind(y,t(fdobjy[[i]]$coefs))}
  # x is the coefficient in the fdobject
  bigx=t(bigDATA)# it is N X (p+1)
  N <- nrow(x)
  p <- ncol(x)
  q<-ncol(y)
  prop <- c()
  # t <- matrix(0,N,K)
  n <- colSums(t)
  prop <- n / N  # props is pi as a vector
  mu <- matrix(NA, K, p)
  for (i in 1:K) mu[i, ] <- colSums(x*t[, i])/n[i]

  ind <- apply(t>0, 2, which)
  n_bis <- c()
  for(i in 1:K) n_bis[i] <- length(ind[[i]])
  #calculation of the degrees of freedom

  #
  #Calculation on Var/Covar matrices for x
  #
  # 
  traceVect <- c()

  # we keep track of the trace (== sum of eigenvalues) to compute the b
  ev <- matrix(0, K, p)
  Q <- vector(mode='list', length=K)
  fpcaobj = list()
  for (i in 1:K){
    donnees <- .mypca.fd2(fdobj, Wlist, t[,i])
    traceVect[i] = sum(diag(donnees$valeurs_propres))
    ev[i, ] <- donnees$valeurs_propres
    Q[[i]] <- donnees$U
    fpcaobj[[i]] = donnees
  }
  # 
  Q1 <- Q
  #Intrinsic dimensions selection
  # browser()
  if (model%in%c("AJBQD", "ABQD")){
    d <- rep(com_dim, length=K)
  } else if ( model%in%c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD") ){
    dmax <- min(apply((ev>noise.ctrl)*rep(1:ncol(ev), each=K), 1, which.max))-1
    if(com_dim>dmax) com_dim <- max(dmax, 1)
    d <- rep(com_dim, length=K)
  } else {
    # 
    d <- .reg_hdclassif_dim_choice(ev, n, method, threshold, FALSE, noise.ctrl)
  }

  #Setup of the Qi matrices

  for(i in 1:K) Q[[i]] <- matrix(Q[[i]][, 1:d[i]], p, d[i])


  #Calculation of the remaining parameters of the selected model

  # PARAMETER a
  ai <- matrix(NA, K, max(d))
  if ( model%in%c('AKJBKQKDK', 'AKJBQKDK', 'AKJBKQKD', 'AKJBQKD') ){
    for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
  } else if ( model%in%c('AKBKQKDK', 'AKBQKDK' , 'AKBKQKD', 'AKBQKD') ){
    for (i in 1:K) ai[i, ] <- rep(sum(ev[i, 1:d[i]])/d[i], length=max(d))
  } else if(model=="AJBQD"){
    for (i in 1:K) ai[i, ] <- ev[1:d[1]]
  } else if(model=="ABQD") {
    ai[] <- sum(ev[1:d[1]])/d[1]
  } else {
    a <- 0
    eps <- sum(prop*d)
    for (i in 1:K) a <- a + sum(ev[i, 1:d[i]])*prop[i]
    ai <- matrix(a/eps, K, max(d))
  }

  # PARAMETER b
  bi <- c()
  denom = min(N,p)
  if ( model%in%c('AKJBKQKDK', 'AKBKQKDK', 'ABKQKDK', 'AKJBKQKD', 'AKBKQKD', 'ABKQKD') ){
    for(i in 1:K){
      remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
      # bi[i] <- sum(ev[i, (d[i]+1):min(N, p)])/(p-d[i])
      bi[i] <- remainEV/(p-d[i]) #pour moi c'est p au lieu de denom
    }
  } else if ( model%in%c("ABQD", "AJBQD") ){
    remainEV = traceVect - sum(ev[1:d[1]])
    # bi[1:K] <- sum(ev[(d[1]+1):min(N, p)])/(min(N, p)-d[1])
    bi[1:K] <- remainEV/(denom-d[1])
  } else {
    b <- 0
    eps <- sum(prop*d)
    for(i in 1:K){
      remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
      # b <- b + sum(ev[i, (d[i]+1):min(N, p)])*prop[i]
      b <- b + remainEV*prop[i]
    }
    bi[1:K] <- b/(min(N, p)-eps)
  }
###########################parameters for regression and covariance for Y
  #gam <- vector(mode='list', length=K)
  #covy<- vector(mode='list', length=K)
  #icovy<- vector(mode='list', length=K)
  ldetcov<- c(rep(1,K))

  ###@@@@@@@@@@@@@@@@@@@@@@to be added as parameters with these default values
  gami=matrix(0, nrow=K, ncol=(p+1)*q )
  covyi    = matrix(0, nrow=K, ncol=q^2 )
  icovyi = matrix(0, nrow=K, ncol=q^2 )
  pqp=p+1
  cres<-.C(".C_mstep",as.character(modely), as.integer(N),as.integer(pqp), as.integer(q),
             as.integer(K), as.double(prop),as.double(bigx),as.double(y),as.double(t), gami = gami,covyi=covyi,
  icovyi=icovyi,logi = as.double(rep(0,K)), as.double(cmtol), as.integer(cmmax), PACKAGE = "funclustweight")
  gam    = array(cres$gami, dim= c(q,p+1,K))
  covy = array(cres$covyi, dim= c(q,q,K) )
  icovy = array(cres$icovyi, dim= c(q,q,K) )
  ldetcov=cres$logi

  

  #---------------------------------------

  list(model = model,
       modely = modely,
       K = K,
       d = d,
       a = ai,
       b = bi,
       mu = mu,
       prop = prop,
       ev = ev,
       Q = Q,
       fpcaobj = fpcaobj,
       Q1 = Q1,
       gam=gam,
       covy=covy,
       icovy=icovy,
       ldetcov=ldetcov
       )


} # end of funclustweight_m_step
#############################################
##################################################

# */

.mypca.fd2 <- function(fdobj,Wlist, Ti){
  if (inherits(fdobj, 'list')){ # multivariate CASE
    #saving the mean before contering
    mean_fd <- list()
    for (i in 1:length(fdobj)){
      mean_fd[[i]] <- fdobj[[i]]
    }
    mean_fd<-list()
    for (i in 1:length(fdobj)){
      mean_fd[[i]]<-fdobj[[i]]
    }

    #centrage des objets fonctionnels
    for (i in 1:length(fdobj)){
      coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj[[i]]$coefs))) * fdobj[[i]]$coefs, 1, sum) / sum(Ti)
      fdobj[[i]]$coefs <- sweep(fdobj[[i]]$coefs, 1, coefmean)
      mean_fd[[i]]$coefs = as.matrix(data.frame(mean=coefmean))
    }



    #Constructing the matrix of coefficents
    coef <- t(fdobj[[1]]$coefs)
    for ( i in 2:length(fdobj) ){
      coef <- cbind( coef, t(fdobj[[i]]$coefs) )
    }
    #covariance matrix

    #Matrice de covariance

    mat_cov <- crossprod( t( .reg_repmat(sqrt(Ti), n = dim(t(coef))[[1]],p=1) *
                               t(coef) ) ) / sum(Ti)
    cov <- Wlist$W_m %*% mat_cov %*% t(Wlist$W_m)
    #eigenvalues and eigenfunctions
    valeurs <- Eigen(cov)
    valeurs_propres <- valeurs$values
    vecteurs_propres <- valeurs$vectors

    bj <- solve(Wlist$W_m) %*% vecteurs_propres
    fonctionspropres <- fdobj[[1]]
    fonctionspropres$coefs <- bj
    scores <- coef %*% Wlist$W %*% bj

    varprop <- valeurs_propres / sum(valeurs_propres)

    pcafd <- list(valeurs_propres = valeurs_propres,
                  harmonic = fonctionspropres,
                  scores = scores,
                  covariance = cov,
                  U = bj,
                  varprop = varprop,
                  meanfd = mean_fd
                  )

  }else if (!inherits(fdobj, 'list')) { #univariate CASE
    #Calculation of means by group
    mean_fd<-fdobj
    #Centrer les objets fonctionnels par groupe
    coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj$coefs))) * fdobj$coefs, 1, sum) / sum(Ti)
    fdobj$coefs <- sweep(fdobj$coefs, 1, coefmean)
    mean_fd$coefs = as.matrix(data.frame(mean=coefmean))


    coef<-t(fdobj$coefs)


    #covariance matrix
    mat_cov <- crossprod( t( .reg_repmat(sqrt(Ti),
                                     n = dim( t(coef) )[[1]],p=1) * t(coef) ) ) / sum(Ti)
    cov <- Wlist$W_m %*% mat_cov %*% t(Wlist$W_m) #A@5 this is the last formula on page 7
    #eigenvalues and eigenfunctions
    #--------------------------------------------
    valeurs <- Eigen(cov)
    valeurs_propres <- valeurs$values
    vecteurs_propres <- valeurs$vectors
    fonctionspropres <- fdobj
    bj <- solve(Wlist$W_m) %*% vecteurs_propres
    fonctionspropres$coefs <- bj
    #calculations of scores by the formula for pca.fd
    scores <- inprod(fdobj, fonctionspropres)

    varprop <- valeurs_propres / sum(valeurs_propres)

    pcafd <- list(valeurs_propres = valeurs_propres,
                  harmonic = fonctionspropres,
                  scores = scores,
                  covariance = cov,
                  U = bj,
                  meanfd = mean_fd
                  )
    #-------------------------------------------
  }
  class(pcafd) <- "pca.fd"
  return(pcafd)
} # end of .mypca.fd2


.reg_hddc_ari <- function(x, y){
  #This function is drawn from the mclust package
  x <- as.vector(x)
  y <- as.vector(y)
  tab <- table(x, y)
  if ( all( dim(tab) == c(1, 1) ) ) return(1)
  a <- sum( choose(tab, 2) )
  b <- sum( choose(rowSums(tab), 2) ) - a
  c <- sum( choose(colSums(tab), 2) ) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
} # end of .reg_hddc_ari



####
#### CONTROLS ####
####

.reg_hddc_control = function(call){

  prefix = "funclustweight: "
  .reg_myCallAlerts(call, "datax", "list,fd", 3, TRUE, prefix)
  .reg_myCallAlerts(call, "datay", "list,fd", 3, TRUE, prefix)
  .reg_myCallAlerts(call, "K", "integerVectorGE1", 3, FALSE, prefix) #
  .reg_myCallAlerts(call, "known", "intNAVectorGE1", 3, FALSE, prefix) #
  .reg_myCallAlerts(call, "model", "characterVector", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "modely", "characterVector", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "threshold", "numericVectorGE0LE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "criterion", "character", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "com_dim", "singleIntegerGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "itermax", "singleIntegerGE2", 3, FALSE, prefix) # IAIN set to 2 iteration min
  .reg_myCallAlerts(call, "eps", "singleNumericGE0", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "graph", "singleLogical", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "d_select", "singleCharacter", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "init", "singleCharacter", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "show", "singleLogical", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "mini.nb", "integerVectorGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "min.individuals", "singleIntegerGE2", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "noise.ctrl", "singleNumericGE0", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "mc.cores", "singleIntegerGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "nb.rep", "singleIntegerGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "keepAllRes", "singleLogical", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "d_max", "singleIntegerGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "d_range", "integerVectorGE1", 3, FALSE, prefix)
  .reg_myCallAlerts(call, "verbose", "singleLogical", 3, FALSE, prefix)


  ####
  #### SPECIFIC controls
  ####

  # Getting some elements
  datax = eval.parent(call[["datax"]], 2)
  datay = eval.parent(call[["datay"]], 2)
  K = eval.parent(call[["K"]], 2)
  known = eval.parent(call[["known"]], 2)
  init = eval.parent(call[["init"]], 2)
  criterion = eval.parent(call[["criterion"]], 2)
  d_select = eval.parent(call[["d_select"]], 2)
  data_length = 0
  N = 0
  d_range = eval.parent(call[["d_range"]], 2)

  if (is.fd(datax)){
    # No NA in the data:
    if (any(is.na(datax$coefs))) stop("funclustweight: NA values in the data are not supported. Please remove them beforehand.", call. = FALSE)
    if (any(!is.numeric(datax$coefs))) stop("funclustweight: Only numeric values are supported in fdata object. Please check coefficients.", call. = FALSE)
    if (any(!is.finite(datax$coefs))) stop("funclustweight: Only finite values are supported in fdata object. Please check coefficients.", call. = FALSE)
    # Size of the data
    if(any(K>2*NROW(datax$coefs))) stop("funclustweight: The number of observations must be at least twice the number of clusters.", call. = FALSE)

    data_length <- nrow(datax$coefs)
    N <- ncol(datax$coefs)

  }else{
    # No NA in the data:
    for (i in 1:length(data)){
      if(!is.fd(datax[[i]])) stop("funclustweight: All dimensions of data must be fd object.", call. = FALSE)
      if (any(is.na(datax[[i]]$coefs))) stop("funclustweight: NA values in the data are not supported. Please remove them beforehand.", call. = FALSE)
      if (any(!is.numeric(datax[[i]]$coefs))) stop("funclustweight: Only numeric values are supported in fdata object. Please check coefficients.", call. = FALSE)
      if (any(!is.finite(datax[[i]]$coefs))) stop("funclustweight: Only finite values are supported in fdata object. Please check coefficients.", call. = FALSE)
      data_length <- data_length + nrow(datax[[i]]$coefs)
    }
    # Size of the data
    if(any(K>2*NROW(datax[[1]]$coefs))) stop("funclustweight: The number of observations must be at least twice the number of clusters ", call. = FALSE)
    N <- ncol(datax[[1]]$coefs)
  }

  if (is.fd(datay)){
    # No NA in the data:
    if (any(is.na(datay$coefs))) stop("funclustweight: NA values in the data are not supported. Please remove them beforehand.", call. = FALSE)
    if (any(!is.numeric(datay$coefs))) stop("funclustweight: Only numeric values are supported in fdata object. Please check coefficients.", call. = FALSE)
    if (any(!is.finite(datay$coefs))) stop("funclustweight: Only finite values are supported in fdata object. Please check coefficients.", call. = FALSE)
    # Size of the data
    if(any(K>2*NROW(datay$coefs))) stop("funclustweight: The number of observations must be at least twice the number of clusters.", call. = FALSE)

    # data_length <- nrow(datay$coefs)
    if(ncol(datay$coefs) != N) stop("funclustweight: Number of curves in X and Y do not match!")

  }else{
    # No NA in the data:
    for (i in 1:length(data)){
      if(!is.fd(datay[[i]])) stop("funclustweight: All dimensions of data must be fd object.", call. = FALSE)
      if (any(is.na(datay[[i]]$coefs))) stop("funclustweight: NA values in the data are not supported. Please remove them beforehand.", call. = FALSE)
      if (any(!is.numeric(datay[[i]]$coefs))) stop("funclustweight: Only numeric values are supported in fdata object. Please check coefficients.", call. = FALSE)
      if (any(!is.finite(datay[[i]]$coefs))) stop("funclustweight: Only finite values are supported in fdata object. Please check coefficients.", call. = FALSE)
      # data_length <- data_length + nrow(datay[[i]]$coefs)
    }
    # Size of the data
    if(any(K>2*NROW(datay[[1]]$coefs))) stop("funclustweight: The number of observations must be at least twice the number of clusters ", call. = FALSE)
    if(ncol(datay[[1]]$coefs) != N) stop("funclustweight: Number of curves in X and Y do not match!")
  }

  # Initialization Controls
  if(!is.null(init)){

    # we get the value of the initialization
    init = .reg_myAlerts(init, "init", "singleCharacterMatch.arg", "funclustweight: ", c('random', 'kmeans', 'tkmeans', 'mini-em', 'vector'))

    # Custom initialization => controls and setup
    if(init == "vector"){
      fdobj = data
      if (!inherits(fdobj, 'list')) {x = t(fdobj$coefs)}
      else {x = t(fdobj[[1]]$coefs); for (i in 2:length(fdobj)) x = cbind(x,t(fdobj[[i]]$coefs))}
      .reg_myCallAlerts(call, "init.vector", "(integer,factor)Vector", 3, FALSE, prefix)

      init.vector = eval.parent(call[["init.vector"]], 2)

      if(is.null(init.vector)) stop("funclustweight: When init='vector', the argument 'init.vector' should be provided.", call. = FALSE)

      if(length(unique(K))>1) stop("funclustweight: Several number of classes K cannot be estimated when init='vector'.", call. = FALSE)

      init.vector <- unclass(init.vector)
      if(K!=max(init.vector)) stop("funclustweight: The number of class K, and the number of classes in the initialization vector are different", call. = FALSE)

      if( length(init.vector)!=nrow(x) ) stop("funclustweight: The size of the initialization vector is different of the size of the data", call. = FALSE)
    }

    # The param init A@5 is this even implemented????
    if (init=='param' && nrow(data)<ncol(data)){
      stop("funclustweight: The 'param' initialization can't be done when N<p", call. = FALSE)
    }

    # The mini.em init
    if (init=='mini-em'){

      mini.nb = eval.parent(call[["mini.nb"]], 2)

      if(!is.null(mini.nb) && length(mini.nb)!=2){
        stop("funclustweight: The parameter mini.nb must be a vector of length 2 with integers\n", call. = FALSE)
      }

    }
  }
  if(!is.null(d_select) && d_select == 'grid') {
    if(!is.null(d_range) && (max(d_range) > data_length)) stop("funclustweight: Intrinsic dimension 'd' cannot be larger than number of input parameters. Please set a lower max.", call. = FALSE)
  }

  if(!is.null(known)) {
    if(all(is.na(known))) stop("funclustweight: declared 'known' must have known values from each class (not all NA).", call. = F)
    if(length(known) != N) stop("funclustweight: 'known' length must match the number of observations from data (known may include NA for observations where groups are unknown).", call. = F)
    if(length(K) > 1) stop("funclustweight: only one group count can be used since known must have values for each group.", call. = F)
    if(length(unique(known[!is.na(known)])) > K) stop("funclustweight: at most K different classes may be passed to 'known'.", call. = F)
    if(max(known, na.rm = T) > K) stop("funclustweight: group numbers must come from integers up to K (ie. for K = 3 integers are from 1, 2, 3).", call. = F)
  }
}

.reg_default_kmeans_control = function(control){

  .reg_myAlerts(control,"kmeans.control","list","kmeans controls: ")

  #
  # Default values of the control parameters
  #

  myDefault = list()
  myDefault$iter.max = 10
  myDefault$nstart = 1
  myDefault$algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen")
  myDefault$trace = FALSE
  myDefault$alpha = 0.2

  #
  # Types of each arg
  #

  myTypes = c("singleIntegerGE1", "singleIntegerGE1", "match.arg",
              "singleLogical", "singleIntegerGE1")

  #
  # Recreation of the kmeans controls + Alerts
  #

  control = .reg_matchTypeAndSetDefault(control, myDefault,
                                    myTypes, "kmeans list of controls: ")

  return(control)
}

#=================================#
# This file contains all the
# "control" functions
#=================================#

# Possible elements of myAlerts:
#
# /* BEGIN FROM FUNHDDC (UNMODIFIED) */
# /*
# * Authors: Bouveyron, C. Jacques, J.
# * Date Taken: 2022-01-01
# * Original Source: funHDDC (unmodified)
# * Address: https://github.com/cran/funHDDC
# *
# */

.reg_myCallAlerts = function(call, name, myType, nParents=1, mustBeThere=FALSE, prefix=""){
  # This function basically calls the function myAlerts, but the arguments are different

  if( name %in% names(call) ){
    # we check the element exists => to provide a fine error
    what = call[[name]]
    val = try(eval.parent(what, nParents), silent = TRUE)
    if( "try-error" %in% class(val) ){
      if( inherits(what, 'name') ){
        # it means the variable was not found
        stop(prefix,"For argument '",name,"': object '",what,"' not found.", call. = FALSE)
      } else {
        stop(prefix,"For argument '",name,"': expression ",as.character(as.expression(what))," could not be evaluated.", call. = FALSE)
      }

    } else {
      a = .reg_myAlerts(val, name, myType, prefix)
      return(a)
    }
  } else if(mustBeThere) {
    stop(prefix, "The argument '", name, "' must be provided.", call. = FALSE)
  }
}

.reg_myAlerts = function(x, name, myType, prefix="", charVec){
  # Format of my types:
  #   - single => must be of lenght one
  #   - Vector => must be a vector
  #   - Matrix => must be a matrix
  #   - GE/GT/LE/LT: greater/lower than a given value
  #   - predefinedType => eg: numeric, integer, etc
  #   - match.arg => very specific => should match the charVec
  # If there is a parenthesis => the class must be of specified types:
  # ex: "(list, data.frame)" must be a list of a data.frame

  ignore.case = TRUE

  firstMsg = paste0(prefix,"The argument '",name,"' ")

  # simple function to extract a pattern
  # ex: if my type is VectorIntegerGE1 => myExtract("GE[[:digit:]]+","VectorIntegerGE1") => 1
  myExtract = function(expr, text, trim=2){
    start = gregexpr(expr,text)[[1]] + trim
    length = attr(start,"match.length") - trim
    res = substr(text,start,start+length-1)
    as.numeric(res)
  }

  #
  # General types handling
  #

  loType = tolower(myType)

  if(grepl("single",loType)){
    if(length(x)!=1) stop(firstMsg,"must be of length one.", call. = FALSE)
  }

  if(grepl("vector",loType) && !grepl("factor",loType)){
    if(!is.vector(x)) stop(firstMsg,"must be a vector.", call. = FALSE)
    if(is.list(x)) stop(firstMsg,"must be a vector (and not a list).", call. = FALSE)
  }

  res = .reg_checkTheTypes(loType, x)
  if(!res$OK) stop(firstMsg,res$message, call. = FALSE)

    # GE: greater or equal // GT: greater than // LE: lower or equal // LT: lower than
  if(grepl("ge[[:digit:]]+",loType)){
    n = myExtract("ge[[:digit:]]+", loType)
    if( !all(x>=n, na.rm = T) ) stop(firstMsg,"must be greater than, or equal to, ", n,".", call. = FALSE)
  }
  if(grepl("gt[[:digit:]]+",loType)){
    n = myExtract("gt[[:digit:]]+", loType)
    if( !all(x>n, na.rm = T) ) stop(firstMsg,"must be strictly greater than ", n,".", call. = FALSE)
  }
  if(grepl("le[[:digit:]]+",loType)){
    n = myExtract("le[[:digit:]]+", loType)
    if( !all(x<=n, na.rm = T) ) stop(firstMsg,"must be lower than, or equal to, ",n,".", call. = FALSE)
  }
  if(grepl("lt[[:digit:]]+",loType)){
    n = myExtract("lt[[:digit:]]+", loType)
    if( !all(x<n, na.rm = T) ) stop(firstMsg,"must be strictly lower than ", n,".", call. = FALSE)
  }

  #
  # Specific Types Handling
  #

  if(grepl("match.arg",loType)){
    if(ignore.case){
      x = toupper(x)
      newCharVec = toupper(charVec)
    } else {
      newCharVec = charVec
    }

    if( is.na(pmatch(x, newCharVec)) ){
      n = length(charVec)
      if(n == 1){
        msg = paste0("'",charVec,"'")
      } else {
        msg = paste0("'", paste0(charVec[1:(n-1)], collapse="', '"), "' or '",charVec[n],"'")
      }
      stop(firstMsg, "must be one of:\n", msg, ".", call. = FALSE)
    } else {
      qui = pmatch(x, newCharVec)
      return(charVec[qui])
    }
  }
}

.reg_matchTypeAndSetDefault = function(myList, myDefault, myTypes, prefix){
  # This  function:
  #   i) check that all the elements of the list are valid
  #   ii) put the default values if some values are missing
  #   iii) Gives error messages if some types are wrong
  # This function obliges  myList to be valid (as given by myDefault)

  # 1) check that the names of the list are valid
  if(is.null(myList)) myList = list()
  list_names = names(myList)

  if(length(list_names)!=length(myList) || any(list_names=="")){
    stop(prefix,"The elements of the list should be named.", call. = FALSE)
  }

  obj_names = names(myDefault)

  isHere = pmatch(list_names,obj_names)

  if(anyNA(isHere)){
    if(sum(is.na(isHere))==1) stop(prefix, "The following argument is not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
    else stop(prefix, "The following arguments are not defined: ",paste(list_names[is.na(isHere)],sep=", "), call. = FALSE)
  }

  # 2) We set the default values and run Warnings
  res = list()
  for(i in 1:length(obj_names)){
    obj = obj_names[i]
    qui = which(isHere==i) # qui vaut le numero de l'objet dans myList
    type = myTypes[i] # we extract the type => to control for "match.arg" type
    if(length(qui)==0){
      # we set to the default if it's missing
      if(type == "match.arg") {
        res[[obj]] = myDefault[[i]][1]
      } else {
        res[[obj]] = myDefault[[i]]
      }
    } else {
      # we check the integrity of the value
      val = myList[[qui]]
      if(type == "match.arg"){
        # If the value is to be a match.arg => we use our controls and not
        # directly the one of the function match.arg()
        charVec = myDefault[[i]]
        .reg_myAlerts(val, obj, "singleCharacterMatch.arg", prefix, charVec)
        val = match.arg(val, charVec)
      } else {
        .reg_myAlerts(val, obj, type, prefix)
      }

      res[[obj]] = val
    }
  }

  return(res)
}



.reg_checkTheTypes = function(str, x){
  # This function takes in a character string describing the types of the
  # element x => it can be of several types

  # types that are controlled for:
  allTypes = c("numeric", "integer", "character", "logical", "list", "data.frame", "matrix", "factor", "intna")

  OK = FALSE
  message = c()

  for(type in allTypes){

    if(grepl(type, str)){
      # we add the type of the control
      if(type == "intna") {
        message = c(message, "integer or NA")
      } else {
        message = c(message, type)
      }

      if(type == "numeric"){
        if(!OK & is.numeric(x)){
          OK = TRUE
        }
      } else if(type == "integer"){
        if(is.numeric(x) && (is.integer(x) || (all(is.finite(x)) && all(x%%1==0)))){ # IAIN added non finite condition
          OK = TRUE
        }
      } else if(type == "intna"){
        if(is.numeric(x[!is.na(x)]) && (is.integer(x[!is.na(x)]) || (all(is.finite(x[!is.na(x)])) && all(x[!is.na(x)]%%1==0)))){ # IAIN added non finite condition
          OK = TRUE
        }
      } else if(type == "character"){
        if(is.character(x)){
          OK = TRUE
        }
      } else if(type == "logical"){
        if(is.logical(x)){
          OK = TRUE
        }
      } else if(type == "list"){
        if(is.list(x)){
          OK = TRUE
        }
      } else if(type == "data.frame"){
        if(is.data.frame(x)){
          OK=TRUE
        }
      } else if(type == "matrix"){
        if(is.matrix(x)){
          OK = TRUE
        }
      } else if(type == "factor"){
        if(is.factor(x)){
          OK = TRUE
        }
      }

    }

    if(OK) break
  }

  if(length(message) == 0) OK = TRUE #ie there is no type to be searched
  else if(length(message) >= 3){
    n = length(message)
    message = paste0("must be of type: ",  paste0(message[1:(n-1)], collapse = ", "), " or ", message[n], ".")
  } else {
    message = paste0("must be of type: ",  paste0(message, collapse = " or "), ".")
  }


  return(list(OK=OK, message=message))
}




.reg_hdclassif_dim_choice <- function(ev, n, method, threshold, graph, noise.ctrl, d_set){
  # Selection of the intrinsic dimension
  # browser()
  N <- sum(n)
  prop <- n/N
  K = ifelse(is.matrix(ev), nrow(ev), 1)

  # browser()

  if(is.matrix(ev) && K>1){
    p <- ncol(ev)
    if(method=="cattell"){
      dev <- abs(apply(ev, 1, diff))
      max_dev <- apply(dev, 2, max, na.rm=TRUE)
      dev <- dev/rep(max_dev, each=p-1)
      d <- apply((dev>threshold)*(1:(p-1))*t(ev[, -1]>noise.ctrl), 2, which.max)

      if(graph){
        # return user settings to original
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        op = par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9), 1+floor(K/4)-1*(K==12)+1*(K==7)))
        for(i in 1:K){
          sub1 <- paste("Class #", i, ",  d", i, "=", d[i], sep="")
          Nmax <- max(which(ev[i, ]>noise.ctrl))-1
          plot(dev[1:(min(d[i]+5, Nmax)), i], type="l", col="blue", main=paste("Cattell's Scree-Test\n", sub1, sep=""), ylab=paste("threshold =", threshold), xlab="Dimension", ylim=c(0, 1.05))
          abline(h=threshold, lty=3)
          points(d[i], dev[d[i], i], col='red')
        }
        par(op)
      }
    } else if(method=="bic"){

      d <- rep(0, K)
      if(graph) op = par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9), 1*(1+floor(K/4)-1*(K==12)+1*(K==7))))

      for (i in 1:K) {
        B <- c()
        Nmax <- max(which(ev[i, ]>noise.ctrl))-1
        p2 <- sum(!is.na(ev[i, ]))
        Bmax <- -Inf
        for (kdim in 1:Nmax){
          if ((d[i]!=0 & kdim>d[i]+10)) break
          a <- sum(ev[i, 1:kdim])/kdim
          b <- sum(ev[i, (kdim+1):p2])/(p2-kdim)
          if (b<0 | a<0){
            B[kdim] <- -Inf
          } else {
            L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
            B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
          }

          if ( B[kdim]>Bmax ){
            Bmax <- B[kdim]
            d[i] <- kdim
          }
        }

        if(graph){
          plot(B, type='l', col=4, main=paste("class #", i, ",  d=", d[i], sep=''), ylab='BIC', xlab="Dimension")
          points(d[i], B[d[i]], col=2)
        }
      }
      if(graph) par(op)
    } else if(method=="grid"){
      d <- d_set
    }
  } else{
    ev <- as.vector(ev)
    p <- length(ev)

    if(method=="cattell"){
      dvp <- abs(diff(ev))
      Nmax <- max(which(ev>noise.ctrl))-1
      if (p==2) d <- 1
      else d <- max(which(dvp[1:Nmax]>=threshold*max(dvp[1:Nmax])))
      diff_max <- max(dvp[1:Nmax])

      if(graph){
        plot(dvp[1:(min(d+5, p-1))]/diff_max, type="l", col="blue", main=paste("Cattell's Scree-Test\nd=", d, sep=''), ylab=paste("threshold =", threshold, sep=' '), xlab='Dimension', ylim=c(0, 1.05))
        abline(h=threshold, lty=3)
        points(d, dvp[d]/diff_max, col='red')
      }
    } else if(method=="bic"){
      d <- 0
      Nmax <- max(which(ev>noise.ctrl))-1
      B <- c()
      Bmax <- -Inf
      for (kdim in 1:Nmax){
        if (d!=0 && kdim>d+10) break
        a <- sum(ev[1:kdim])/kdim
        b <- sum(ev[(kdim+1):p])/(p-kdim)
        if (b<=0 | a<=0) B[kdim] <- -Inf
        else{
          L2 <- -1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))*N
          B[kdim] <- 2*L2 - (p+kdim*(p-(kdim+1)/2)+1)*log(N)
        }
        if ( B[kdim]>Bmax ){
          Bmax <- B[kdim]
          d <- kdim
        }
      }

      if(graph){
        plot(B, type='l', col=4, main=paste("BIC criterion\nd=", d, sep=''), ylab='BIC', xlab="Dimension")
        points(d, B[d], col=2)
      }
    }
  }
  return(d)
}
################################
.Reg_hdclassift_bic <- function(par, p,q){
  model <- par$model
  modely <-par$modely
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  N <- par$N
  prop <- par$prop

  if(length(b)==1){
    #update of b to set it as variable dimension models
    eps <- sum(prop*d)
    n_max <- if(model%in%c("ABQD", "AJBQD")) length(par$ev) else ncol(par$ev)
    b <- b*(n_max-eps)/(p-eps)
    b <- rep(b, length=K)
  }
  if (length(a)==1) a <- matrix(a, K, max(d))
  else if (length(a)==K) a <- matrix(a, K, max(d))


  if(min(a, na.rm=TRUE)<=0 | any(b<0)) return(-Inf)

  if (is.null(par$loglik)){
    som_a <- c()
    for (i in 1:K) som_a[i] <- sum(log(a[i, 1:d[i]]))
    L <- -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop) + (p+q)*(1+log(2*pi))))*N
  }  else L <- par$loglik[length(par$loglik)]


  ro <- K*p+K-1
  tot <- sum(d*(p-(d+1)/2))
  D <- sum(d)
  d <- d[1]
  to <- d*(p-(d+1)/2)
  if (model=='AKJBKQKDK') m <- ro+tot+D+K
  else {if (model=='AKBKQKDK') m <- ro+tot+2*K
  else if (model=='ABKQKDK') m <- ro+tot+K+1
  else if (model=='AKJBQKDK') m <- ro+tot+D+1
  else if (model=='AKBQKDK') m <- ro+tot+K+1
  else if (model=='ABQKDK') m <- ro+tot+2}
  ###############
  if (modely=='EII') m <- m+1
  else {if (modely=='VII') m <- m+K
  else if (modely=='EEI') m <- m+q
  else if (modely=='VEI') m <- m+K+q-1
  else if (modely=='EVI') m <- m+1+K*(q-1)
  else if (modely=='VVI') m <- m+K*q
  else if (modely=='EEE') m <- m+q*(q+1)/2
  else if (modely=='VEE') m <- m+K+q-1+q*(q-1)/2
  else if (modely=='EVE') m <- m+1+K*(q-1)+q*(q-1)/2
  else if (modely=='EEV') m <- m+q+K*q*(q-1)/2
  else if (modely=='VVE') m <- m+q*K+q*(q-1)/2
  else if (modely=='VEV') m <- m+K+q-1+K*q*(q-1)/2
  else if (modely=='EVV') m <- m+1+K*(q-1)+K*q*(q-1)/2
  else if (modely=='VVV') m <- m+K*q*(q+1)/2}
  bic <- -(-2*L+m*log(N))

  #calcul ICL
  t = par$posterior
  #if(!is.null(t)){
  # means we are in HDDC
  Z = ((t - apply(t, 1, max))==0) + 0
  icl = bic - 2*sum(Z*log(t+1e-15))
  # } else {
  #   # Si HDDA, entropie est nulle => car classes pures
  #   icl = bic
  #}

  return(list(bic = bic, icl = icl))
}

#########################################
.Reg.hdc_getComplexity = function(par, p,q){
  model <- par$model
  modely <- par$modely
  K <- par$K
  d <- par$d
  b <- par$b
  a <- par$a
  mu <- par$mu
  N <- par$N
  prop <- par$prop

  ro <- K*p+K-1
  tot <- sum(d*(p-(d+1)/2))
  D <- sum(d)
  d <- d[1]
  to <- d*(p-(d+1)/2)
  if (model=='AKJBKQKDK') m <- ro+tot+D+K
  else {if (model=='AKBKQKDK') m <- ro+tot+2*K
  else if (model=='ABKQKDK') m <- ro+tot+K+1
  else if (model=='AKJBQKDK') m <- ro+tot+D+1
  else if (model=='AKBQKDK') m <- ro+tot+K+1
  else if (model=='ABQKDK') m <- ro+tot+2}
  if (modely=='EII') m <- m+1
  else {if (modely=='VII') m <- m+K
  else if (modely=='EEI') m <- m+q
  else if (modely=='VEI') m <- m+K+q-1
  else if (modely=='EVI') m <- m+1+K*(q-1)
  else if (modely=='VVI') m <- m+K*q
  else if (modely=='EEE') m <- m+q*(q+1)/2
  else if (modely=='VEE') m <- m+K+q-1+q*(q-1)/2
  else if (modely=='EVE') m <- m+1+K*(q-1)+q*(q-1)/2
  else if (modely=='EEV') m <- m+q+K*q*(q-1)/2
  else if (modely=='VVE') m <- m+q*K+q*(q-1)/2
  else if (modely=='VEV') m <- m+K+q-1+K*q*(q-1)/2
  else if (modely=='EVV') m <- m+1+K*(q-1)+K*q*(q-1)/2
  else if (modely=='VVV') m <- m+K*q*(q+1)/2}

  return(m)
}

######################
.reg_hdc_getTheModely = function(model, all2models = FALSE){
  # Function used to get the models from number or names

  model_in = model

  if(!is.vector(model)) stop("The argument 'model' must be a vector.")

  if(anyNA(model)) stop("The argument 'model' must not contain any NA.")

  ModelNames <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV","EVV","VVV")

  model = toupper(model)

  if(length(model)==1 && model=="ALL"){
    if(all2models) model <- 1:14
    else return("ALL")
  }

  qui = which(model %in% 1:14)
  model[qui] = ModelNames[as.numeric(model[qui])]

  # We order the models properly
  qui = which(!model%in%ModelNames)
  if (length(qui)>0){
    if(length(qui)==1){
      msg = paste0("(e.g. ", model_in[qui], " is incorrect.)")
    } else {
      msg = paste0("(e.g. ", paste0(model_in[qui[1:2]], collapse=", or "), " are incorrect.)")
    }
    stop("Invalid model name ", msg)
  }

  # warning:
  if(max(table(model))>1) warning("The model vector, argument 'model', is made unique (repeated values are not tolerated).")

  mod_num <- c()
  for(i in 1:length(model)) mod_num[i] <- which(model[i]==ModelNames)
  mod_num <- sort(unique(mod_num))
  model <- ModelNames[mod_num]

  return(model)
}


##############################

.reg_hdc_getTheModel = function(model, all2models = FALSE){
  # Function used to get the models from number or names

  model_in = model

  if(!is.vector(model)) stop("The argument 'model' must be a vector.")

  if(anyNA(model)) stop("The argument 'model' must not contain any NA.")

  ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK")

  model = toupper(model)

  if(length(model)==1 && model=="ALL"){
    if(all2models) model <- 1:6
    else return("ALL")
  }

  qui = which(model %in% 1:6)
  model[qui] = ModelNames[as.numeric(model[qui])]

  # We order the models properly
  qui = which(!model%in%ModelNames)
  if (length(qui)>0){
    if(length(qui)==1){
      msg = paste0("(e.g. ", model_in[qui], " is incorrect.)")
    } else {
      msg = paste0("(e.g. ", paste0(model_in[qui[1:2]], collapse=", or "), " are incorrect.)")
    }
    stop("Invalid model name ", msg)
  }

  # warning:
  if(max(table(model))>1) warning("The model vector, argument 'model', is made unique (repeated values are not tolerated).")

  mod_num <- c()
  for(i in 1:length(model)) mod_num[i] <- which(model[i]==ModelNames)
  mod_num <- sort(unique(mod_num))
  model <- ModelNames[mod_num]

  return(model)
}



####
#### Utilities ####
####


.reg_addCommas = function(x) sapply(x, .reg_addCommas_single )

.reg_addCommas_single = function(x){
  # This function adds commas for clarity for very long values of likelihood

  if(!is.finite(x)) return(as.character(x))

  s = sign(x)
  x = abs(x)

  decimal = x - floor(x)
  if(decimal>0) dec_string = substr(decimal, 2, 4)
  else dec_string = ""

  entier = as.character(floor(x))

  quoi = rev(strsplit(entier, "")[[1]])
  n = length(quoi)
  sol = c()
  for(i in 1:n){
    sol = c(sol, quoi[i])
    if(i%%3 == 0 && i!=n) sol = c(sol, ",")
  }

  res = paste0(ifelse(s==-1, "-", ""), paste0(rev(sol), collapse=""), dec_string)
  res
}


.reg_repmat <- function(v,n,p){ #A@5 WHAT IS THIS???????????????????????
  if (p==1){M = cbind(rep(1,n)) %*% v} #A@5 a matrix of column of v
  else { M = matrix(rep(v,n),n,(length(v)*p),byrow=T)} # removed cat("!");
  M
}

.reg_diago <- function(v){
  if (length(v)==1){ res = v }
  else { res = diag(v)}
  res
}

# /* END OF FROM FUNHDDC */
#####################################
C_mahalanobis<-function(bigx, y, gami, icovyi){
  regsub <-t(y)-gami%*%t(bigx)
  # IAIN suggest to transpose-ok
  inter<-t(regsub)%*%icovyi%*%regsub
  auxint <-diag(inter)
  return (auxint)
}

#########################################C version
.reg_imahalanobis <- function(x, muk, wk, Qk, aki) {
  # w = fpcaobj[[i]]$W
  # *************************************************************************** #
  # should be called as imahalanobis(x, mu[i,], w[i,], Q[[i]], a[i,], b[i])
  # return formula on top of page 9
  # *************************************************************************** #

  # Code modified to take vectors instead of matrices
  p <- ncol(x)
  N <- nrow(x)
  res <- rep(0, N)



  #X <- x - matrix(muk, N, p, byrow=TRUE)

  #  Qi <- wk %*% Qk


  # xQi <- (X %*% Qi)

  #proj <- (X %*% Qi) %*% aki
  #the R result- no C code
  # res_old <- rowSums(proj ^ 2)



  # Calling C_imahalanobis
  res <- .C(".C_imahalanobis", as.double(x), as.double(muk), as.double(wk),
            as.double(Qk), as.double(aki), as.integer(p),
            as.integer(N), as.integer(nrow(aki)), res = rep(0,N), PACKAGE="funclustweight")


  return(res$res)

} # end of .reg_imahalanobis

#********************************* TIMING FUNCTION *********************************#
# /*
# * Authors: Andrews, J. Wickins, J. Boers, N. McNicholas, P.
# * Date Taken: 2023-01-01
# * Original Source: teigen (modified)
# * Address: https://github.com/cran/teigen
# *
# */
.reg_estimateTime <- function(stage, start.time, totmod, backspace = 0){
  curwidth <- ifelse(backspace != 0, backspace, getOption("width"))
  ##Output enough backspace characters to erase the previous output progress information
  ##\b is used instead of \r due to RStudio's difficulty with handling carriage returns
  # cat(paste(rep("\b", backspace), collapse = ""))

  ##the string that will eventually be output to the screen
  output.string <- ""

  ##[short,med,long]string all are strings that will be output depending on size of the
  ##console. These strings will be one of two things, depending on if stage is "init" or not
  if(stage=="init"){
    medstring <- longstring <- "????"
    shortstring <- "0"
    modelcount <- NA
  }
  else{
    for(i in 1:(backspace)) cat("\b")
    # cat("\f")

    modelcount <- stage+1 ##number of models calculated
    modsleft <- totmod-modelcount
    timerun <- difftime(Sys.time(),start.time, units="secs")
    timeremain <- (timerun/modelcount)*modsleft

    if(timeremain>60){
      units(timeremain) <- "mins"
    } else if(timeremain>3600){
      units(timeremain) <- "hours "
    } else if(timeremain>86400){
      units(timeremain) <- "days"
    } else if(timeremain>604800){
      units(timeremain) <- "weeks "
    }

    if(timerun>60){
      units(timerun) <- "mins"
    } else if(timerun>3600){
      units(timerun) <- "hours "
    } else if(timerun>86400){
      units(timerun) <- "days"
    } else if(timerun>604800){
      units(timerun) <- "weeks "
    }

    ## % complete
    shortstring <- paste(format(round((1-modsleft/totmod)*100), width=3, justify="right"), sep = "")
    ##approx. time remaining
    medstring <- paste(format(round(timeremain,1), width=4, justify="right"), sep = "")
    ##Time taken
    longstring <- paste(format(round(timerun,1), width=4, justify="right"), sep = "")
  }

  ##start to build up the output string depending on console size
  if(curwidth>=15){
    shortstring <- str_pad(shortstring, 5, pad=" ")
    output.string <- paste(shortstring, "% complete", sep="")

    ##if larger console, output more information
    if(curwidth>=48){
      medstring <-str_pad(medstring, 10, pad=" ")
      output.string <- paste("Approx. remaining:", medstring,"  |  ", output.string, sep = "")

      ##if even larger console, output even more information
      if(curwidth>=74){
        longstring <- str_pad(longstring, 10, pad=" ")
        output.string <- paste("Time taken:", longstring, "  |  ", output.string, sep = "")
      }
    }
  }

  cat(output.string)
  flush.console()
  ## return model count and output string size
  c(modelcount, nchar(output.string))


} # end of .reg_estimateTime

######################################
###########################################
