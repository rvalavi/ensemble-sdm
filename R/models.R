# install required packages
pkg <- c(
    "dismo",
    "gbm",
    "mgcv",
    "glmnet",
    "terra",
    "blockCV",
    "dismo",
    "precrec"
)

pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
if (length(pkgna)) {
    nm <- paste(pkgna, collapse = ", ")
    message("Installing packages: ", nm, "...")
    utils::install.packages(pkgna)
}


# remotes::install_github("rvalavi/myspatial")
# library(myspatial)
suppressPackageStartupMessages({
    library(dismo)
    library(gbm)
    library(mgcv)
    library(glmnet)
    library(dismo) #***** or predicts?
    library(precrec)
})


# change the formula of the variables with low unique value (this causes issue in GAM)
update_names <- function(x, y, k = 4, thr = 20) {
    z <- as.numeric(
        which(
            apply(x[, y, drop = FALSE], 2, function(t) length(unique(t)) < thr)
        )
    )
    
    for (i in z) {
        y[i] <- paste0(y[i], ", k = ", k)
    }
    
    return(y)
}


# Fitting several models as an ensemble SDM
#
# A function to ... 
# Models include: GAM, Lasso-GLM, RF, BRT (a.k.a. GBM), and Maxent
#
# x data.frame of the training data.
# y character, the column name of the response
# use_ranger logical, use ranger instead of randomForest for the RF component
# balance_classes logical, draw equal numbers of both classes per RF tree
ensemble <- function(
        x,
        y = "",
        fold_ids = NULL,
        models = c("GLM", "GAM", "GBM", "RF", "Maxent"),
        use_ranger = FALSE,
        balance_classes = TRUE
) {
    
    x <- as.data.frame(x)
    # get covariate names
    covars <- names(x)
    covars <- covars[covars != y]
    
    if(is.null(fold_ids)){
        fold_ids <- dismo::kfold(x[, y, drop = TRUE], k = 5)
    }
    
    # if null use all
    if (is.null(models)) models <- c("GLM", "GAM", "GBM", "RF", "Maxent")
    if (identical(tolower(models)[1], "all")) models <- c("GLM", "GAM", "GBM", "RF", "Maxent")
    models <- tolower(models)
    if (!any(models %in% c("glm", "gam", "gbm", "rf", "maxent"))) {
        stop("Select at least one supported ensemble model.", call. = FALSE)
    }
    
    # set the default to null
    ls_mod <- gm_mod <- br_mod <- rf_mod <- mx_mod <- quad_obj <- NULL
    
    # compute weights
    num_pr <- as.numeric(table(x[, y, drop=TRUE])["1"]) # number of presences
    num_bg <- as.numeric(table(x[, y, drop=TRUE])["0"]) # number of backgrounds
    case_weights <- ifelse(x[, y, drop = TRUE] == 1, 1, num_pr / num_bg)
    
    if ("glm" %in% models) {
        message("Fitting GLM-Lasso...")
        # fit regularized regression with L1
        quad_obj <- make_quadratic(x, cols = covars)
        training_quad <- predict(quad_obj, newdata = x)
        new_vars <- names(training_quad)[names(training_quad) != y]
        training_sparse <- Matrix::sparse.model.matrix(~., training_quad[, new_vars, drop = FALSE])
        
        ls_mod <- glmnet::cv.glmnet(
            x = training_sparse,
            y = training_quad[, y],
            family = "binomial",
            alpha = 1,
            weights = case_weights,
            foldid = fold_ids
        )
        plot(ls_mod)
    }
    
    
    if ("gam" %in% models) {
        require(mgcv)
        message("Fitting GAM...")
        # fit generalised additive model (GAM))
        
        # update name with k for case that have low unique values and cause
        # issue in smoothing
        var_names <- update_names(x, covars, k = 4, thr = 20)
        # building the model formula here to include updated names
        gm_form <- reformulate(
            termlabels = paste0("s(", var_names, ")"), 
            response = y
        )
        
        gm_mod <- mgcv::gam(
            formula = as.formula(gm_form),
            data = x,
            family = binomial(link = "logit"),
            weights = case_weights,
            method = "REML",
            select = TRUE
        )
    }
    
    
    if ("gbm" %in% models) {
        for (learning_rate in c(0.001, 0.0005, 0.00025)) {
            message("Fitting GBM (learning rate ", learning_rate, ")...")
            br_mod <- dismo::gbm.step(
                data = x,
                gbm.y = which(names(x) == y),
                gbm.x = which(names(x) != y),
                family = "bernoulli",
                tree.complexity = ifelse(num_pr < 50, 1, 5),
                learning.rate = learning_rate,
                bag.fraction = 0.75,
                max.trees = 10000,
                n.trees = 50,
                site.weights = case_weights,
                fold.vector = fold_ids,
                n.folds = length(unique(fold_ids)),
                silent = TRUE
            )
            if (!is.null(br_mod)) break
        }
        if (is.null(br_mod)) {
            stop("GBM fitting failed at all three learning rates.", call. = FALSE)
        }
    }
    
    
    # for RF mtry
    sqrt_vars <- floor(
        sqrt(
            length(covars)
        )
    )
    # RF for binary classification
    if ("rf" %in% models) {
        rf_mod <- random_forest(
            data = x,
            y = y,
            splitrule = c("gini", "hellinger"),
            num.trees = 1000,
            mtry = seq.int(min(2, length(covars)), min(length(covars), max(5, sqrt_vars))),
            foldid = fold_ids,
            balance_classes = balance_classes,
            use_ranger = use_ranger
        )
    }
    
    
    # NOTE: Maxent models is only for binary presence-background data
    if ("maxent" %in% models) {
        message("Fitting Maxent...")
        # tune maxent parameters
        param_optim <- maxent_param(
            data = x,
            y = y,
            folds = fold_ids
        )
        print(param_optim)
        
        # fit a maxent model with the tuned parameters
        mx_mod <- dismo::maxent(
            x = x[, covars, drop = FALSE],
            p = x[[y]],
            args = param_optim
        )
        
    }
    
    out <- list(
        "GLM-Lasso" = ls_mod,
        "GAM" = gm_mod,
        "GBM" = br_mod,
        "RF" = rf_mod,
        "Maxent" = mx_mod,
        "quad" = quad_obj
    )
    class(out) <- c("ensemble", "list")

    # Fit the scale on training predictions and reuse it for every projection.
    training_predictions <- predict(out, x[, covars, drop = FALSE], .components = TRUE)
    valid <- colSums(!is.finite(training_predictions)) == 0L
    if (any(!valid)) {
        stop(
            "Non-finite training predictions from: ",
            paste(colnames(training_predictions)[!valid], collapse = ", "),
            call. = FALSE
        )
    }
    attr(out, "predictors") <- covars
    attr(out, "prediction_ranges") <- apply(training_predictions, 2, range)

    return(out)
}

# print method for ensemble object
print.ensemble <- function(x, ...){
    print(class(x))
    mods <- names(x)[!sapply(x, is.null)]
    mods <- mods[mods != "quad"]
    cat("Models:", paste(mods, collapse = ", "), "\n")
}

# predict with ensemble object
predict.ensemble <- function(object, newdata, type = "response", ..., .components = FALSE){
    newdata <- as.data.frame(newdata)
    ranges <- attr(object, "prediction_ranges")
    if (!is.null(ranges)) {
        if (!identical(type, "response")) {
            stop("The scaled ensemble only supports type = 'response'.", call. = FALSE)
        }
        newdata <- newdata[, attr(object, "predictors"), drop = FALSE]
        complete <- complete.cases(newdata)
        result <- rep(NA_real_, nrow(newdata))
        if (!any(complete)) return(result)
        newdata <- newdata[complete, , drop = FALSE]
    }

    pred_ls <- pred_gm <- pred_rf <- pred_br <- pred_mx <- NULL

    # predict lasso with the minimum-error lambda
    if (!is.null(object[["GLM-Lasso"]])) {
        pred_quad <- predict(object[["quad"]], newdata = newdata)
        pred_sparse <- Matrix::sparse.model.matrix(~., pred_quad)
        pred_ls <- predict(object[["GLM-Lasso"]], pred_sparse, s = "lambda.min", type = type, ...)[, 1]
    }

    # predict gam model
    if (!is.null(object[["GAM"]])) {
        pred_gm <- predict(object[["GAM"]], newdata, type = type, ...)
    }

    # predict brt with the best tree number
    brt <- object[["GBM"]]
    if (!is.null(brt)) {
        pred_br <- predict(brt, newdata, n.trees = brt$gbm.call$best.trees, type = type, ...)
    }

    # predict rf with the fitted backend
    rf <- object[["RF"]]
    if (inherits(rf, "ranger")) {
        loadNamespace("ranger")
        pred_rf <- predict(rf, newdata, type = type, num.threads = 5, ...)$predictions[, "1"]
    } else if (!is.null(rf)) {
        loadNamespace("randomForest")
        pred_rf <- predict(rf, newdata, type = "prob", ...)[, "1"]
    }

    # predict maxent on the cloglog scale
    mxt <- object[["Maxent"]]
    if (!is.null(mxt)) {
        pred_mx <- dismo::predict(mxt, newdata, args = "outputformat=cloglog")
    }

    out <- cbind(
        "GLM-Lasso" = pred_ls,
        "GAM" = pred_gm,
        "GBM" = pred_br,
        "RF" = pred_rf,
        "Maxent" = pred_mx
    )

    if (.components) return(out)
    if (is.null(ranges)) return(rowMeans(out, na.rm = TRUE))

    out[!is.finite(out)] <- NA_real_
    spans <- ranges[2, ] - ranges[1, ]
    constant <- spans == 0
    spans[constant] <- 1
    scaled <- sweep(sweep(out, 2, ranges[1, ], "-"), 2, spans, "/")
    # Use the midpoint when a component has a constant training score.
    scaled[, constant] <- ifelse(is.na(scaled[, constant, drop = FALSE]), NA_real_, 0.5)
    result[complete] <- rowMeans(pmax(pmin(scaled, 1), 0))
    return(result)
}


# function for simultaneous tuning of maxent regularization multiplier and features
maxent_param <- function(data, y = "occ", folds = NULL, k = 5, filepath = tempdir()){
    if (!isTRUE(dismo::maxent(silent = TRUE))) {
        stop("MaxEnt requires rJava and a working Java runtime.", call. = FALSE)
    }

    if(is.null(folds)){
        # generate balanced CV folds
        folds <- dismo::kfold(data[, y, drop = TRUE], k = k)
    }
    covars <- names(data)[which(names(data) != y)]
    names(data)[which(names(data) == y)] <- "occ"
    # regularization multipliers
    ms <- c(0.5, 1, 2, 3, 4)
    grid <- expand.grid(
        regmult = paste0("betamultiplier=", ms),
        features = list(
            c("noautofeature", "nothreshold"), # LQHP
            c("noautofeature", "nothreshold", "noproduct"), # LQH
            c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
            c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
            c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
        stringsAsFactors = FALSE
    )
    AUCs <- rep(NA_real_, nrow(grid))
    nfold <- sort(unique(folds))
    last_error <- NULL
    for (n in seq_len(nrow(grid))) {
        parameters <- as.character(unlist(grid[n, ]))
        scores <- rep(NA_real_, nrow(data))

        for (k in nfold) {
            train_set <- which(folds != k)
            test_set <- which(folds == k)
            fit_path <- tempfile("maxent-", tmpdir = filepath)
            fold_scores <- tryCatch({
                maxmod <- dismo::maxent(
                    x = data[train_set, covars, drop = FALSE],
                    p = data$occ[train_set],
                    path = fit_path,
                    args = parameters
                )
                dismo::predict(
                    maxmod, data[test_set, covars, drop = FALSE], args = "outputformat=cloglog"
                )
            }, error = function(error) error, finally = unlink(fit_path, recursive = TRUE))
            if (inherits(fold_scores, "error")) {
                last_error <- conditionMessage(fold_scores)
                break
            }
            scores[test_set] <- fold_scores
        }

        if (all(is.finite(scores))) {
            AUCs[n] <- precrec::auc(precrec::evalmod(scores = scores, labels = data$occ))[1, 4]
        }
    }
    if (all(is.na(AUCs))) {
        detail <- if (is.null(last_error)) "" else paste0(" Last error: ", last_error)
        stop("No MaxEnt candidate completed cross-validation.", detail, call. = FALSE)
    }
    best_param <- as.character(unlist(grid[which.max(AUCs), ]))

    return(best_param)
}


# Tune and refit the selected random forest backend.
random_forest <- function(data, 
                          y = "occ", 
                          splitrule = c("gini", "hellinger", "extratrees"),
                          num.trees = 500,
                          mtry = NULL,
                          # max.depth = NULL,
                          foldid = NULL,
                          probability = TRUE,
                          sample.fraction = 0.632,
                          case.weights = NULL, 
                          type = "response",
                          threads = NULL,
                          plot = TRUE,
                          balance_classes = FALSE,
                          use_ranger = TRUE){
    if (!is.logical(use_ranger) || length(use_ranger) != 1L || is.na(use_ranger)) {
        stop("use_ranger must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(balance_classes) || length(balance_classes) != 1L || is.na(balance_classes)) {
        stop("balance_classes must be TRUE or FALSE.", call. = FALSE)
    }

    package <- if (use_ranger) "ranger" else "randomForest"
    if (!length(find.package(package, quiet = TRUE))) {
        message("Installing packages: ", package, "...")
        utils::install.packages(package)
    }
    message("Fitting RF with ", package, "...")

    predictors <- setdiff(names(data), y)
    form <- as.formula(
        paste(y, "~ .")
    )
    if (probability) {
        data[, y] <- as.factor(data[, y])
    }
    
    
    if(is.null(mtry)){
        mtry <- floor(sqrt(ncol(data) - 1))
    }
    
    grid <- expand.grid(split = if (use_ranger) splitrule else "gini",
                        ntrees = num.trees,
                        mtry = mtry,
                        # depth = max.depth, 
                        stringsAsFactors = FALSE)
    
    if(is.null(foldid)){
        foldid <- dismo::kfold(data[, y, drop = TRUE], k = 5)
    }
    nfold <- sort(unique(foldid))
    
    fit_forest <- function(training, params, weights) {
        response <- training[[y]]
        if (balance_classes) {
            presences <- sum(response == "1")
            backgrounds <- sum(response == "0")
            if (presences == 0L || backgrounds == 0L) {
                stop("Balanced RF fitting requires presences and backgrounds in every training fold.",
                     call. = FALSE)
            }
        }

        if (!use_ranger) {
            if (balance_classes) {
                # randomForest caps each stratum's sample size at its observed count.
                draws <- min(presences, backgrounds)
                sampsize <- c("0" = draws, "1" = draws)
            } else {
                sampsize <- max(1L, floor(sample.fraction * nrow(training)))
            }
            return(randomForest::randomForest(
                x = training[, predictors, drop = FALSE],
                y = response,
                ntree = params$ntrees,
                mtry = params$mtry,
                sampsize = sampsize,
                weights = weights,
                replace = TRUE
            ))
        }

        # ranger uses fractions of the entire training set for each class.
        fraction <- if (balance_classes) rep(presences / nrow(training), 2) else sample.fraction
        ranger::ranger(
            formula = form,
            data = training,
            num.trees = params$ntrees,
            mtry = params$mtry,
            splitrule = params$split,
            probability = probability,
            sample.fraction = fraction,
            case.weights = weights,
            num.threads = threads,
            replace = TRUE
        )
    }

    evalmodel <- grid
    evalmodel$eval <- NA_real_
    for(i in seq_along(grid[,1])){
        modauc <- c()
        for(k in nfold){
            
            train_set <- which(foldid != k)
            test_set <- which(foldid == k)
            
            mod <- fit_forest(data[train_set, , drop = FALSE], grid[i, ], case.weights[train_set])
            test_data <- data[test_set, predictors, drop = FALSE]
            if (use_ranger) {
                pred <- predict(mod, test_data, type = type)$predictions[, "1"]
            } else {
                pred <- predict(mod, test_data, type = "prob")[, "1"]
            }
            modauc <- c(modauc, precrec::auc(
                precrec::evalmod(scores = pred, labels = data[test_set, y, drop = TRUE])
            )[1, 4])
        }
        evalmodel$eval[i] <- mean(modauc)
    }
    
    bestparam <- which.max(evalmodel$eval)
    
    # if(plot){
    #     print(
    #         plot_tune_param(evalmodel)
    #     )
    # }
    print(evalmodel[bestparam, ])
    
    final_model <- fit_forest(data, grid[bestparam, ], case.weights)

    return(final_model)
}

# plot tuning parameters of ranger
plot_tune_param <- function(x){
    levs <- sort(unique(x$depth))
    x$depth <- as.character(x$depth)
    x$depth <- ifelse(x$depth == 0, "max", x$depth)
    x$depth <- factor(x$depth, levels = c(as.character(levs[-1]), "max"))
    x$id <- as.numeric(as.factor(x$mtry))
    return(
        ggplot(data = x, aes(x = id, y = eval, col = depth, shape = split)) +
            geom_point(aes(size = as.factor(mtry))) +
            geom_path() +
            facet_grid(~ntrees) +
            scale_x_continuous(breaks = unique(x$id)) +
            theme_bw() +
            labs(x = "", y = "R-squared", size = "mtry")
    )
}



# Orthogonal quadratic polynomials for glmnet
#
# A function to create quadratic terms for glmnet functions i.e. lasso and ridge regression.
# The output is an object of make_quadratic that can be used to predict on rasters and data.frames
# for creating the quadratic terms.
make_quadratic <- function(df, cols = NULL, verbose = TRUE){
    if(is.null(cols)){
        cols <- colnames(df)
    }
    if(is.numeric(cols)){
        cols <- colnames(df)[cols]
    }
    if(!all(is.element(cols, colnames(df)))){
        stop("The cols should be the same as the column names.")
    }
    # remove the factors
    factors <- vapply(df[cols], is.factor, logical(1))
    if(any(factors) && verbose){
        message("The factor columns were removed from cols: ", paste(cols[factors], collapse = ", "))
    }
    cols <- cols[!factors]
    xbar <- vapply(df[cols], mean, numeric(1))
    x1 <- data.frame(mapply(`-`, df[cols], xbar, SIMPLIFY = FALSE))
    sum_squares <- colSums(x1 ^ 2)
    sum_squares[sum_squares == 0] <- 1
    alpha <- colSums(x1 ^ 3) / sum_squares
    # specify the output class
    finalList <- list(names = cols, xbars = xbar, alphas = alpha)
    class(finalList) <- c("make_quadratic")
    return(finalList)
}


predict.make_quadratic <- function(object, newdata, ...){
    if(!methods::is(object, "make_quadratic"))
        stop("object should be a make_quadratic object.")
    if(!all(object$names %in% names(newdata)))
        stop("The newdata does not have the same names as the object.")
    ncl <- object$names
    if(methods::is(newdata, "Raster")){
        for(i in ncl){
            x1 <- newdata[[i]] - object$xbars[i]
            x2 <- (x1 ^ 2) - (object$alphas[i] * x1)
            if(raster::nlayers(newdata) > 1){
                newdata <- newdata[[-which(names(newdata) == i)]]
                newdata <- raster::stack(newdata, x1)
            } else{
                newdata <- x1
            }
            names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_1")
            newdata <- raster::stack(newdata, x2)
            names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_2")
        }
    } else if(methods::is(newdata, "data.frame")){
        for(i in ncl){
            x1 <- newdata[,i] - object$xbars[i]
            x2 <- x1 ^ 2 - object$alphas[i] * x1
            newdata <- newdata[,-which(names(newdata) == i), drop = FALSE]
            newdata[,ncol(newdata) + 1] <- x1
            names(newdata)[ncol(newdata)] <- paste0(i, "_1")
            newdata[,ncol(newdata) + 1] <- x2
            names(newdata)[ncol(newdata)] <- paste0(i, "_2")
        }
    } else stop("newdata should be a raster or a data.frame.")
    return(newdata)
}

