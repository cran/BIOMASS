correctTaxo = function (genus, species = NULL, score = 0.5) 
{
  options(stringsAsFactors = F)
  
  # Create a dataframe with the original values
  oriData <- data.frame(genus = as.character(genus), query = as.character(genus), id = 1:length(genus))
    
  # Create a temporary vector with only the unique of the names
  query <- na.omit(unique(genus))
  
  if(!is.null(species)) 
  {
    # Check the length of the inputs
    if(length(genus) != length(species))
      stop("You should provide two vectors of genus and species of the same length")
    
    # Create a dataframe with the original values
    oriData <- data.frame(genus = as.character(genus), species = as.character(species), 
                          query = paste(genus, species), id = 1:length(genus))
    
    # Create a temporary vector with only the unique of the names
    query <- unique(paste(genus, species))
    query <- na.omit(query)
  }
  
  if (length(query) < 1 || is.na(query)) 
    stop("Please supply at least one name", call. = FALSE)

  getpost <- "get"
  if(length(query) > 50)
    getpost <- "post"
  
  tc <- function(l) Filter(Negate(is.null), l)
  con_utf8 <- function(x) httr::content(x, "text", encoding = "UTF-8")
  
  # If there is too much data, better submit it in separated queries
  splitby <- 30
  slicedQu <- rep(1, length(query))
  if (getpost == "get" && length(query) > 75 | length(query) > 30 && getpost == "post") 
    slicedQu <- rep(1:ceiling(length(query) / splitby), each = splitby)[1:length(query)]
  
  url <- "http://taxosaurus.org/submit"
  
  # Send and retrieve the data from taxosaurus
  df <- c()
  for(s in unique(slicedQu))
  {
    x <- query[slicedQu == s]
    
    if(getpost == "get") 
    {
      query2 <- paste(gsub(" ", "+", x, fixed = T), collapse = "%0A")
      args <- tc(list(query = query2))
      out <- httr::GET(url, query = args)
      retrieve <- out$url
    }
    else 
    {
      loc <- tempfile(fileext = ".txt")
      write.table(data.frame(x), file = loc, col.names = FALSE, row.names = FALSE)
      args <- tc(list(file = httr::upload_file(loc), source = "iPlant_TNRS"))
      out <- httr::POST(url, body = args, httr::config(followlocation = 0))
      tt <- con_utf8(out)
      message <- jsonlite::fromJSON(tt, FALSE)[["message"]]
      retrieve <- jsonlite::fromJSON(tt, FALSE)[["uri"]]
    }
    
    print(paste("Calling", retrieve))
    
    timeout <- "wait"
    while (timeout == "wait") 
    {
      ss <- httr::GET(retrieve)
      output <- jsonlite::fromJSON(con_utf8(ss), FALSE)
      if(!grepl("is still being processed", output["message"]) == TRUE) 
        timeout <- "done"
    }
    out <- tc(output$names)
    
    for(j in 1:length(out))
    {
      if(length(out[[j]][[2]][[1]]$annotations) == 0)
        out[[j]][[2]][[1]]$annotations <- list(Authority = NA)
      
      tmp <- data.frame(out[[j]][[2]][[1]])
      tmp$submittedName <- out[[j]]$submittedName
      df <- rbind(df, tmp)
    }
  }
  
  # Remove some parasite characters 
  df$submittedName <- gsub("\"", "", df$submittedName)
  df$submittedName <- gsub("\r", "", df$submittedName)
  df$matchedName <- gsub("\"", "", df$matchedName)
  df$matchedName <- gsub("\r", "", df$matchedName)
  
  #### AUTOMATIC PROCEDURE
  # If score ok
  df$outName[df$score >= score] <- df$matchedName[df$score >= score]
  
  # If score non ok
  df$outName[df$score < score] <- df$submittedName[df$score < score]
    
  df <- unique(df[, c("submittedName", "outName")])
  df <- df[with(df, order(submittedName)), ]
  
  df$nameModified <- TRUE
  df$nameModified[!is.na(df$outName) & df$outName == df$submittedName] <- FALSE
  
  out <- merge(oriData, df, by.x = "query", by.y = "submittedName", all.x = T)
  out$nameModified[is.na(out$nameModified)] <- TRUE
  out <- out[with(out, order(id)), ]
  
  out$genusCorrected <- sapply(strsplit(out$outName, "[ ]"), "[", 1)
  out$speciesCorrected <- sapply(strsplit(out$outName, "[ ]"), "[", 2)
  
  return(out[, c("genusCorrected", "speciesCorrected", "nameModified")])
}