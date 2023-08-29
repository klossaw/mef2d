

liftCellChat_ADAPTED <- function(object, group.new = NULL) {
  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    if (is.null(group.new)) {
      group.max.all <- unique(unlist(sapply(idents, levels)))
      group.num <- sapply(idents, nlevels)
      group.num.max <- max(group.num)
      group.max <- levels(idents[[which(group.num == group.num.max)]])
      if (length(group.max) != length(group.max.all)) {
        stop("CellChat object cannot lift up due to the missing cell groups in any dataset. Please define the parameter `group.new`!")
      }
    } else {
      group.max <- group.new
      group.num.max <- length(group.new)
    }
    message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    for (i in 1:length(idents)) {
      cat("Update slots object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      group.i <- levels(idents[[i]])
      group.existing <- group.max[group.max %in% group.i]
      # CHANGED THIS LINE: in order to get the indices of the cell types of group.new that exist in the orignal levels we changed it to this: 
      group.existing.index <- unlist(lapply(group.i, function(x) which(group.max %in% x)))
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                              dimnames = list(group.max, group.max))
          values.new[group.existing.index, group.existing.index] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("pairwiseRank")) {
          for (k in 1:length(values)) {
            values.new1 <- vector("list", group.num.max)
            values.new1[group.existing.index] <- values[[k]]
            temp <- values[[k]][[1]]
            temp$prob  <- 0; temp$pval <- 1
            for (kk in setdiff(1:group.num.max, group.existing.index)) {
              values.new1[[kk]] <- temp
            }
            names(values.new1) <- group.max
            values[[k]] <- values.new1
          }
          values.new <- vector("list", group.num.max)
          values.new[group.existing.index] <- values
          temp <- lapply(values.new1, function(x) {
            x$prob  <- 0; x$pval <- 1
            return(x)
          })
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new[[kk]] <- temp
          }
          names(values.new) <- group.max
        }
        net[[net.j]] <- values.new
      }
      object@net[[i]] <- net
      
      # cat("Update slot object@netP...", '\n')
      netP <- object@netP[[i]]
      for (netP.j in names(netP)) {
        values <- netP[[netP.j]]
        if (netP.j %in% c("pathways")) {
          values.new <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("prob")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("centr")) {
          for (k in 1:length(values)) {
            values.new <- lapply(values, function(x) {
              values.new2 <- lapply(x, function(x) {
                values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
                values.new1[group.existing.index] <- x
                names(values.new1) <- group.max
                return(values.new1)
              })
              names(values.new2) <- names(x)
              return(values.new2)
            })
            names(values.new) <- names(values)
          }
        }
        netP[[netP.j]] <- values.new
      }
      object@netP[[i]] <- netP
      # cat("Update slot object@idents...", '\n')
      # idents[[i]] <- factor(group.max, levels = group.max)
      idents[[i]] <- factor(idents[[i]], levels = group.max)
    }
    object@idents[1:(length(object@idents)-1)] <- idents
  } else {
    if (is.null(group.new)) {
      stop("Please define the parameter `group.new`!")
    } else {
      group.max <- as.character(group.new)
      group.num.max <- length(group.new)
      message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    }
    cat("Update slots object@net, object@netP, object@idents in a single dataset...", '\n')
    # cat("Update slot object@net...", '\n')
    net <- object@net
    idents <- object@idents
    group.i <- levels(idents)
    group.existing <- group.max[group.max %in% group.i]
    # CHANGED THIS LINE: in order to get the indices of the cell types of group.new that exist in the orignal levels we changed it to this: 
    group.existing.index <- unlist(lapply(group.i, function(x) which(group.max %in% x)))
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                            dimnames = list(group.max, group.max))
        values.new[group.existing.index, group.existing.index] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("pairwiseRank")) {
        for (k in 1:length(values)) {
          values.new1 <- vector("list", group.num.max)
          values.new1[group.existing.index] <- values[[k]]
          temp <- values[[k]][[1]]
          temp$prob  <- 0; temp$pval <- 1
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new1[[kk]] <- temp
          }
          names(values.new1) <- group.max
          values[[k]] <- values.new1
        }
        values.new <- vector("list", group.num.max)
        values.new[group.existing.index] <- values
        temp <- lapply(values.new1, function(x) {
          x$prob  <- 0; x$pval <- 1
          return(x)
        })
        for (kk in setdiff(1:group.num.max, group.existing.index)) {
          values.new[[kk]] <- temp
        }
        names(values.new) <- group.max
      }
      net[[net.j]] <- values.new
    }
    object@net <- net
    
    
    # cat("Update slot object@netP...", '\n')
    netP <- object@netP
    for (netP.j in names(netP)) {
      values <- netP[[netP.j]]
      if (netP.j %in% c("pathways")) {
        values.new <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("prob")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("centr")) {
        for (k in 1:length(values)) {
          values.new <- lapply(values, function(x) {
            values.new2 <- lapply(x, function(x) {
              values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
              values.new1[group.existing.index] <- x
              names(values.new1) <- group.max
              return(values.new1)
            })
            names(values.new2) <- names(x)
            return(values.new2)
          })
          names(values.new) <- names(values)
        }
      }
      netP[[netP.j]] <- values.new
    }
    object@netP <- netP
    
    # cat("Update slot object@idents...", '\n')
    idents <- factor(idents, levels = group.max)
    object@idents <- idents
  }
  
  return(object)
}


