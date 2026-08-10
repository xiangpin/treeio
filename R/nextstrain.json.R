#' @title read.nextstrain.json
#' @param x the json tree file of auspice from nextstrain.
#' @return treedata object
#' @export
#' @author Shuangbin Xu
#' @examples
#' file1 <- system.file("extdata/nextstrain.json", "minimal_v2.json", package="treeio") 
#' tr <- read.nextstrain.json(file1)
#' tr
read.nextstrain.json <- function(x){
    file1 <- x
    x <- jsonlite::read_json(x)
    if (all(c('meta', 'tree') %in% names(x))){
        dt <- parser_children(x$tree)
    }else{
        dt <- parser_children(x)
    }
    tr <- build_nextstrain_treedata(dt)
    tr@file <- file1
    return(tr)
}

parser_children <- function(x, id=list2env(list(id = 0L)), parent = 1){
    rows <- list()
    stack <- list(list(node = x, parent = parent))
    node_id <- 0L
    while (length(stack) > 0L) {
        cur <- stack[[length(stack)]]
        stack[[length(stack)]] <- NULL
        node_id <- node_id + 1L
        row <- extract_node_attrs(cur$node, id = node_id, parent = cur$parent)
        if ('children' %in% names(cur$node)) {
            kids <- cur$node$children
            if (length(kids) > 0L) {
                for (i in rev(seq_along(kids))) {
                    stack[[length(stack) + 1L]] <- list(node = kids[[i]], parent = node_id)
                }
            }
        }
        rows[[node_id]] <- row
    }
    dat <- dplyr::bind_rows(rows)
    numeric_cols <- vapply(dat, check_num, logical(1))
    if (any(numeric_cols)) {
        dat[numeric_cols] <- lapply(dat[numeric_cols], as.numeric)
    }
    if ('div' %in% colnames(dat)){
        dat[["branch.length"]] <- dat[["div"]] - dat[["div"]][match(dat[["parent"]], dat[["node"]])]
    }
    return(dat)
}

check_num <- function(x){
    is_numeric(x) && is.character(x)
}

extract_node_attrs <- function(x, id, parent){
    if ('node_attrs' %in% names(x)){
        res <- build_node_attrs(x[['node_attrs']])
    }else if('attr' %in% names(x)){
        res <- build_node_attrs(x[['attr']])
    }else{
        res <- data.frame()
    }
    if ('name' %in% names(x)){
        res$name <- x[['name']]
    }else if('strain' %in% names(x)){
        res$name <- x[['strain']]
    }
    res$parent <- parent
    res$node <- id
    return(res)
}

build_node_attrs <- function(x){
    x <- unlist(x)
    index <- grepl('\\.value$', names(x))
    names(x)[index] <- gsub('\\.value$', '', names(x)[index])
    x <- as.data.frame(as.list(x), stringsAsFactors = FALSE, check.names = FALSE)
    return(x)
}

#' @importFrom tidytree tip.label node.label tip.label<- node.label<-
build_nextstrain_treedata <- function(x){
    if ("branch.length" %in% colnames(x)){
        clnm <- c('parent', 'node', 'branch.length')
	phylo <- as.phylo(x[, clnm, drop = FALSE], 
                          branch.length =  'branch.length')
    }else{
        clnm <- c('parent', 'node')
        phylo <- as.phylo(x[, clnm, drop = FALSE])
    }
    tmptbl <- as_tibble(phylo)
    tmptbl$label <- as.integer(tmptbl$label) 
    x <- dplyr::left_join(tmptbl[, !colnames(tmptbl) %in% clnm[clnm != 'node']], x, by = c('label' = 'node')) |>
         dplyr::arrange(!!rlang::sym("node"))
    if ("name" %in% colnames(x)){
        tip.label(phylo) <- x$name[seq(Ntip(phylo))]
        tmpnm <- x$name[seq(Ntip(phylo) + 1, Nnode(phylo, internal.only = FALSE))]
        if (all(is.na(tmpnm)) || all(tmpnm == "")){
            node.label(phylo) <- NULL
        }else{
            node.label(phylo) <- tmpnm
        }
        clnm <- c(clnm, c("name", "label"))
    }    
    trda <- new("treedata", phylo = phylo)
    x <- x[, !colnames(x) %in% clnm[clnm != 'node'], drop=FALSE]
    if (ncol(x)>1){
        trda@data <- x
    } 
    return(trda)
}

