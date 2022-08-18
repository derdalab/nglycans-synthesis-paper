# R script to read campaign file, load LiGA data, match sequences and sort
# with dictionary and report and output processed results - EC 2021-09-15


# extract the various description fields from the LiGA datafile names
# inputs: a dataframe and the index of a column of the filenames
# results are additional columns prepended to the result frame
parse_LiGA_filenames <- function(files, index) {
    date_regex <- "(2[0-9]{3})(0[1-9]|1[012])(0[1-9]|[12][0-9]|3[01])"
    chem_regex <- "-([0-9]+)([A-Z]{2,})([a-z]{2,})([A-Z]{2,})-([A-Z]+)"
    parts <- str_match(basename(as.character(files[[index]])),
                       paste0(date_regex, chem_regex))
    parsed_names <- data.frame(parts[,2:9], files, stringsAsFactors = FALSE)
    colnames(parsed_names)[1:8] <- c("Year", "Month", "Day", "Library",
                                 "Modification", "Target", "Postprocess", "ID")
    parsed_names
}


# given the root (or more) of the data filename and a directory:
# find the file's full name, then read it into a data.frame
# handles both space- and tab-separated files and gzip compression
read.LiGA.data <- function(filename, dir = "."){
    # search for the given filename part and obtain the full filename
    # obtain the full list of files and filter with grepl - avoids a regex
    
    filename.list <- list.files(path = dir)
    matches.index <- which(grepl(pattern = trimws(filename),
                                 filename.list, fixed = TRUE))
    if(length(matches.index) == 0){
        stop(paste0("No file matching '", filename, "' found in '", dir, "'."))
    } else if(length(matches.index) > 1){
        warning(paste0("Found multiple files, using '",
                       filename.list[matches.index[1]], "'."))
    }
    full.name <- file.path(dir, filename.list[matches.index[1]])
    print(paste0('Reading data from file: ', full.name))

    # handle case of a zip file - use first matching filename in ZIP directory
    if(grepl("\\.zip$", full.name)){
        zip.directory <- unzip(full.name, list = TRUE)
        matches.index <- which(grepl(pattern = trimws(filename),
                                     zip.directory$Name, fixed = TRUE))
        zip.name <- zip.directory$Name[matches.index[1]]
        zip.length = zip.directory$Length[matches.index[1]]
        # decompress and copy file from ZIP archive to temporary file
        temp.file.name <- tempfile()
        source.connection <- unz(full.name, zip.name, "rb")
        writeBin(readBin(source.connection, "raw", zip.length, 1),
                 temp.file.name, 1, useBytes = TRUE)
        close(source.connection)
        # switch the data source to the temporary file
        full.name <- temp.file.name
    }

    # read in start of the file as a list of lines and find header line
    file.header.lines <- readLines(full.name, n = 100)
    header.line <- which(grepl(" Mod Nuc AA ", file.header.lines))[1]

    # read in a data file - try reading as a tab separated table first
    file.data <- tryCatch(
        read.table(full.name, header = TRUE, skip = header.line - 1,
                   sep = "\t", quote = "\"", stringsAsFactors = FALSE),
        error = function(condition) {NULL} # silence errors if not tab-separated
    )
    
    # if the input fails or read with only one column, retry as space separated
    if(is.null(file.data) || NCOL(file.data) == 1){
        file.data <- tryCatch(
            read.table(full.name, header = TRUE, skip = header.line - 1,
                       sep = " ", quote = "\"", stringsAsFactors = FALSE),
            error = function(condition) {NULL}  # again silence errors
        )
    }

    # if read.table has failed - try the following line parsing instead
    # for space separated tables when the Modification column is unquoted
    # and contains spaces
    if(is.null(file.data)){
        # read in the file as a list of lines
        file.lines <- readLines(full.name)
      
        # find indices of the header line and the Mod column within the header
        header.line <- which(grepl("Nuc", file.lines))[1]
        mod.column <- which(grepl("Mod",
                                strsplit(file.lines[header.line], " ")[[1]]))
        header <- strsplit(file.lines[header.line], " ")[[1]]
        post.mod.columns <- length(header) - mod.column

        # split each line using spaces
        # keep the leftmost mod - 1 columns and the rightmost post.mod.columns
        # merge the remaining column(s) back into one column using spaces,
        # this is the Mod column value
        resplit.lines <- lapply(
            strsplit(file.lines[(header.line + 2):length(file.lines)], " "),
            function(x){
                len <- length(x);
                c(x[1:(mod.column - 1)],
                  paste(x[mod.column:(len-post.mod.columns)], collapse = " "),
                  x[(len-post.mod.columns+1):len]
                 )
            }
        )

	# copy to data.frame, converting any columns of only numbers to numeric
        file.data <- data.frame(row.names = 1:length(resplit.lines))
        for(i in 1:header.columns){
            data.column <- unlist(lapply(resplit.lines, function(x){x[i]}))
            conversion <- suppressWarnings(as.numeric(data.column))
            if(!any(is.na(conversion))){  # if column is all numeric convert it
                data.column <- conversion
            }
            file.data$data <- data.column
            colnames(file.data)[i] <- header[i]
        }
    }
    
    # as sanity check, ensure that there is a Nuc (nucleotide sequence) column
    if(!"Nuc" %in% colnames(file.data)){
         warning(paste0("Expected column Nuc not found.  File '", full.name,
                        "' may not be read correctly"))
    }

    # delete temporary file if it was used
    if(exists("temp.file.name")){
        unlink(temp.file.name)
    }
    
    # return the data.frame
    file.data
}


### read in the dictionary file, keep (non-blank) SDB and Axis.name pairs
### remove leading and trailing periods in the column names,
### and duplicate SDB entries
read_dictionary_tables <- function(dictionary_basename, order_table_basename,
        label_column = "Axis.name", known_SDB_list, dict_dir = "dictionaries",
        order_table_dir = "order-tables") {

    # list possible dictionaries - if more than one is found, use the first
    dict_filenames <- list.files(path = dict_dir,
                        pattern = paste0(dictionary_basename, ".txt"))
    if(length(dict_filenames) < 1){
        stop("Dictionary not found: ", dictionary_basename)
    }

    # read in the choosen file
    dictionary <- read.delim(file.path(dict_dir, dict_filenames[1]),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # clean-up by stripping leading and trailing periods from column names
    colnames(dictionary) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(dictionary))

    # for simplicity drop the columns we are not using
    dictionary <- dictionary[,colnames(dictionary) %in% c("SDB", label_column, "Alphanum")]

    # standardize the name of the label column
    colnames(dictionary)[colnames(dictionary) == label_column] <- "Axis.name"

    # handle duplicated SDB entries - keep only one copy and label it "mix"
    dupl_dictionary_SDBs <- duplicated(dictionary$SDB)
    dictionary <- dictionary[!dupl_dictionary_SDBs,]
    dictionary$Axis.name[dictionary$SDB %in% dupl_dictionary_SDBs] = "mix"

    # drop unknown SDBs and blank entries from the dictionary
    dictionary <- dictionary[dictionary$SDB %in% known_SDB_list,]
    dictionary <- dictionary[dictionary$Axis.name != "" & dictionary$SDB != "",]

    # check that there is still at least one entry and three columns
    if(length(dictionary) != 3 | length(dictionary$SDB) < 1){
       stop(paste0("The dictionary file appears to be empty or lacks the SDB, Alphanum, and/or Label columns: ",
           dict_filenames[1], " ", input_columns$Label[i,]))
    }

    # read in the order table - drop extra columns and duplicate Alphanum entries
    order_table <- read.delim(file.path(order_table_dir, paste0(order_table_basename, ".csv")), 
                              header = TRUE, sep = ",", stringsAsFactors = FALSE)
    colnames(order_table) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(order_table))
    order_table <- order_table[,colnames(order_table) %in% c("Alphanum", "Order")]

    # remove all copies of duplicate Alphanum entries - multiples propagate through joins
    # TODO: could rescue repeated Alphanum's where all entries agree on the same order
    dupl_order_table_Alphanum <- order_table$Alphanum[duplicated(order_table$Alphanum)]
    order_table <- order_table[!order_table$Alphanum %in% dupl_order_table_Alphanum,]

    if(length(order_table) != 2 | length(order_table$Order) < 1){
        stop(paste0("The Order Table file appears to be empty or lacks the Order column: ",
                    order_table_basename))
    }

    # join the order table to the dictionary and reduce it to
    # the SDB, Axis Label, and Order columns
    dictionary <- dictionary %>% left_join(order_table, by = "Alphanum")
    dictionary <- dictionary[,colnames(dictionary) %in% c("SDB", "Axis.name", "Order")]
    
    # make NA in the Order a special position at the end of the list
    # and any duplicated SDBs a special order position at the front of the list
    dictionary$Order[is.na(dictionary$Order)] <- max(dictionary$Order, na.rm = TRUE) + 1000
    dictionary$Order[dictionary$SDB %in% dupl_dictionary_SDBs] <- min(dictionary$Order) - 1

    # sort the dictionary based on the Order column
    dictionary <- dictionary[order(dictionary$Order, nchar(dictionary$SDB), dictionary$SDB),]

    # return the dictionary table
    dictionary
}

### read file into dataframe listing datafiles, columns to select, and matching dictionary
### return a single merged table of the requested columns
read_and_merge_columns <- function(input_columns, LiGA.data.dir = "LiGA-data",
        dict_dir = "dictionaries", order_table_dir = "order-tables") {

    # verify the inputs
    if(length(colnames(input_columns)[colnames(input_columns) %in%
               c("Filename", "Columns", "Dictionary", "Labels", "OrderTable")]) != 5){
        stop("One or more required columns (Filename, Columns, Dictionary, Labels, OrderTable) not found.")
    }
    # TODO check that columns are numeric values

    if(length(unique(input_columns$OrderTable)) != 1){
        stop("Multiple OrderTables in input .")
    }


    # TODO: what if filenames do not meet the expected pattern?
    parsed_names <- parse_LiGA_filenames(input_columns, "Filename")
    parsed_names$root <- unite(unite(unite(parsed_names, Date, Year:Day, sep = ""),
                       root, Library:Postprocess, sep = ""), 
                       filename_root, Date:ID, sep = "-")$filename_root

    # if an ExcludedSDBs column is not present create an empty one
    if(!"ExcludedSDBs" %in% colnames(parsed_names)){
       parsed_names$ExcludedSDBs <- ""
    }
    if(!"Type" %in% colnames(parsed_names)){
       parsed_names$Type <- ""
    }
    if(!"Reference" %in% colnames(parsed_names)){
        parsed_names$Reference <- ifelse(parsed_names$Type == "Test", "Control", "")
    }

    # read in the table of pre-labelled SDB sequences
    identified_seqs <- read.table("./identified-sequence-table-extra.gz",
                                  header = TRUE, stringsAsFactors = FALSE)

    # create data structures to store the results
    merged_data <- data.frame(SDB = character(), Barcode = character(),
        stringsAsFactors = FALSE)
    column_types <- c("Barcode", "Barcode")
    type.list <- unique(as.character(parsed_names$Type))
    data.source <- data.frame(Filename = c("", ""), Column = c(NA, NA), stringsAsFactors = FALSE)
    data.types <- setNames(data.frame(matrix(FALSE, nrow = 2, ncol = length(type.list)), stringsAsFactors = FALSE), type.list)

    # loop over the listed files, merging columns into the output data frame
    for(i in 1:length(parsed_names$Filename)){
        # read in the file
        file_data <- read.LiGA.data(parsed_names$Filename[i], dir = LiGA.data.dir)
###     file_data <- read.table(parsed_names$Filename[i], header = TRUE,
###               skip = 0, sep = "\t", quote = "\"", stringsAsFactors = FALSE)

        # discard sequences that are not 96 base in length
        #file_data <- file_data[nchar(file_data$Nuc) == 96,]
        # discard the sum total line (Nuc == XX) - can be recalculated later
        file_data <- file_data[file_data$Nuc != "XX",]

        # match each sequence with known SDBs
        file_data <- full_join(file_data, identified_seqs, by = "Nuc")
        noncount.columns <- colnames(file_data) %in%
                c("Distance", "index", "mindex", "Primer", "Nuc", "Mod", "AA")
        file_data <- file_data[,!noncount.columns]

        # handle any missing values - zero out in count columns
        # replace any NAs in the count data with zeros - restore NA in SDB column
        sdb_names <- file_data$SDB
        file_data[is.na(file_data)] <- 0
        file_data$SDB <- sdb_names
        file_data <- file_data %>% group_by(SDB) %>% summarize_all(sum)

        # read in the dictionary file and order table
        dictionary <- read_dictionary_tables(parsed_names$Dictionary[i], parsed_names$OrderTable[i],
                          parsed_names$Label[i], known_SDB_list = identified_seqs$SDB,
                          dict_dir = dict_dir, order_table_dir = order_table_dir)

        # use Axis.name labels where the SDBs match, fallback to the SDB if unknown
        # and drop data for any excluded cases
        file_data <- full_join(file_data, dictionary, by = "SDB")
        file_data$Barcode <- file_data$Axis.name
        unknown_barcode <- is.na(file_data$Axis.name)
        file_data$Barcode[unknown_barcode] <- "" #file_data$SDB[unknown_barcode]
        file_data <- file_data[!file_data$SDB %in% parsed_names$ExcludedSDBs[[i]],]
        file_data <- file_data[order(file_data$Order, nchar(file_data$SDB), file_data$SDB, na.last = TRUE),]

        # list the indices of columns containing the read counts
        # exclude columns named "Mod", "Nuc", "AA", etc.
        data_columns <- (1:NCOL(file_data))[!colnames(file_data) %in% c("SDB",
            "Distance", "index", "mindex", "Primer", "Nuc", "Mod", "AA",
            "Axis.name", "Order", "Barcode")]
        # choose the requested subset of the data columns
        active_columns <- data_columns[as.vector(parsed_names$Columns[[i]])]
        for(col in 1:length(active_columns)){
            old.row <- data.source$Filename == parsed_names$Filename[i] &
                       data.source$Column == active_columns[col]
            if(!any(old.row)){
                # if not already read, append the requested column/row to the count data
                data.source <- rbind(data.source, data.frame(
                        Filename = parsed_names$Filename[i], Column = active_columns[col]))
                data.types <- rbind(data.types, data.types[1,])
            }
            # now that the column is added to data.types mark this entry
             data.types[data.source$Filename == parsed_names$Filename[i] &
                       data.source$Column == active_columns[col],
                                                  as.character(parsed_names$Type[i])] <- TRUE
            # read in the column if new
            if(!any(old.row)){
                column_types <- c(column_types,
                                  as.character(parsed_names$Type[i]))
                if(active_columns[col] > ncol(file_data)){
                   stop("Requested column not in data file. Column:",
                        active_columns[col])
                }
                # new column name: filename root followed by the column number
                # column number as in the campaign file
                column_name <- paste0(parsed_names$root[i], 
                                      parsed_names$Columns[[i]][col])
                colnames(file_data)[active_columns[col]] <- column_name

                # extract a data_frame with the barcode and one column
                # of read counts, summing read counts over barcodes
                data_column <- file_data %>% group_by(SDB, Barcode) %>%
                           summarize(count_sum = sum(!!as.name(column_name)),
                                     .groups = 'drop')
                # note: added experimental .groups option to avoid warning

                # fix the column name and merge it to the result frame
                colnames(data_column)[3] <- column_name
                merged_data <- full_join(merged_data, data_column,
                                         by = c("SDB", "Barcode"))

                # replace any NAs in the count data with zeros - restore NA SDBs
                sdb_names <- merged_data$SDB
                merged_data[is.na(merged_data)] <- 0
                merged_data$SDB <- sdb_names
            }
        }
    }
    ordering_frame <- merged_data %>% left_join(dictionary,
                                        by = c("SDB", "Barcode" = "Axis.name"))
    merged_data$Order <- ordering_frame$Order
    data.types <- rbind(data.types, data.types[1,])
    ordering <- order(ordering_frame$Order,		# primary sort by provided Order
                      ordering_frame$Barcode == "",		# empty labels at the end
                      sub("-\\[[0-9<]*]$", "",  ordering_frame$Barcode),   # break ties by glycan (alphabetic)
                      suppressWarnings(as.integer(
                          sub(".*-\\[<?0*([1-9][0-9]*)\\]$", "\\1", ordering_frame$Barcode))),
                      ordering_frame$Barcode,		# finally just use the raw string
                      substr(ordering_frame$SDB,1,1),
                      nchar(ordering_frame$SDB),		# then use the SDB
                      ordering_frame$SDB, na.last = TRUE)
#       colnames(order_table) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(order_table))
#ddd_frame <<- ordering_frame
    merged_data <- merged_data[ordering,]
# print(ordering_frame[ordering, c("Order", "SDB", "Barcode")])

    # finally extract a table of requested two-way comparisons
    contrasts <- unique(parsed_names[parsed_names$Reference != "", c("Type","Reference")])

    list(table = merged_data, Contrasts = contrasts, Column.types = data.types, dict = dictionary)
}

make.barplot <- function(output_df, labels.on, basename, font.family = font.family) {
    
    # filter out the unknown SDBs - EC 2022-07-07
    #####output_df <- output_df[!grepl("^AzOH", output_df$Barcode) & output_df$Barcode != "",]
    ####output_df <- output_df[!grepl("^AzOH", output_df$Barcode) & output_df$Barcode != "",]
    output_df <- output_df[output_df$Barcode != "",]

    # for the plot estimate the standard deviation of the model output FC
    relative.error <- sqrt(output_df$Test_SD ** 2    / output_df$Test_CPM ** 2 +
                           output_df$Control_SD ** 2 / output_df$Control_CPM ** 2)
    plot_data <- data.frame(
                Label = reorder(factor(paste0(output_df$Barcode, " ", output_df$SDB)), output_df$Order),
                Order = output_df$Order,
                Fold.change = 2 ** output_df$logFC,
                QValue = p.adjust(output_df$PValue, "BH"),
                Noise = (2 ** output_df$logFC) * (1 + relative.error),
                Valid = output_df$Valid
                )

    #plot_data <- plot_data[!is.na(plot_data$QValue),]
    error.bar.data <- plot_data[!is.na(plot_data$Noise),]  # avoid plotting NAs
    significant.data.pos <- error.bar.data[error.bar.data$Fold.change > 1 & error.bar.data$QValue <= 0.05 & error.bar.data$Valid,]
    significant.data.neg <- error.bar.data[error.bar.data$Fold.change < 1 & error.bar.data$QValue <= 0.05 & error.bar.data$Valid,]

    bar_plot <- ggplot(data = plot_data,
                     aes(x = Label, y = Noise, fill = Valid)) +
            theme_light() +
            theme(text = element_text(family = font.family),
                  axis.title = element_text(family = font.family, size = 7),
                  legend.position = "top",
                  legend.title = element_text(size = 7),
                  legend.text = element_text(size = 7),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank()) +
            ggtitle(unique(paste0(filename_root, ": ", difference.frame$Type, " vs. ", difference.frame$Reference))) +
            labs(x = "Glycan", y = "Fold Change") +
            #scale_y_sqrt(limits = c(0, 158.99), expand = expand_scale(mult = c(0, 0.0))) +
            #scale_y_continuous(limits = c(0, max(plot_data$Noise)),
            #                   expand = expand_scale(mult = c(0, 0.10))) +
            scale_y_sqrt(expand = expand_scale(mult = c(0, 0.08))) +
            #scale_y_continuous(trans = "identity",
            #                   # limits = NA, #c(0, 15.0),
            #                   expand = expand_scale(mult = c(0, 0.08))) +
            #scale_y_sqrt() +
            #scale_y_log10() +
            geom_hline(yintercept = 1) + 
            geom_col(aes(y = Fold.change), position = "dodge", size = 0.0, width = 0.60) +
            geom_errorbar(data = error.bar.data,	# top of error bars
                          aes(ymin = Noise, ymax = Noise, width = 0.80)) +
            geom_linerange(data = error.bar.data,	# vertical of error bar
                           aes(ymin = Fold.change, ymax = Noise))
    if(nrow(significant.data.pos) > 0){
        bar_plot <- bar_plot + geom_text(label = "*", hjust = 0.5,
                                  data = significant.data.pos, aes(y = Noise))
    }
    
    #if(nrow(significant.data.neg) > 0){
    #    bar_plot <- bar_plot + geom_text(label = "*", hjust = 0.5,
    #                data = significant.data.neg, color = "red", aes(y = Noise))
    #}

    bar_plot <- bar_plot +
        theme(axis.text.y  = element_text(face = "bold", size = 7)) +
        scale_fill_manual(values = c("#aaaaaa", "#000000"),
                          #labels = c("FALSE", "TRUE"),
                          #guide = guide_legend(reverse = TRUE)
                          )
            theme(axis.text.x  = element_text(face = "bold", angle = 90, hjust = 1, vjust = 0.5, size = 7))

    if(labels.on){
        bar_plot <- bar_plot + theme(axis.text.x  = element_text(
                   face = "bold", angle = 90, hjust = 1, vjust = 0.5, size = 7))
    } else {
        bar_plot <- bar_plot + theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank())
    }

    bar_plot = bar_plot
}

