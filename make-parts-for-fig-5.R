# R script to read campaign file, load LiGA data match and sort with
# dictionary and report and output processed results - EC 2021-09-15


campaign.filename.pattern <- "in-vitro.csv"


LiGA.data.dir   <- "LiGA-data"
dict_dir        <- "dictionaries"
order_table_dir <- "order-tables"
campaign_dir    <- "campaigns"
output.dir      <- "outputs"		# TODO - create if not yet existing
heatmap.filename <- "heatmap.pdf"

font.family     <- "Arial"  # likely options "Arial", "Helvetica", "Arimo"
                            # must be installed on system
                            # mainfont in header of report.Rmd should match 

library(dplyr)
library(tidyverse)
library(knitr)
library(rmarkdown)
library(edgeR)
library(scales)
library(grid)
library(gridExtra)
library(gtable)
source("read-inputs.R")
source("colour-scale.R")



if (!dir.exists(output.dir)){
    dir.create(output.dir, recursive = TRUE)
}

# height and width are in inches - width choosen for a two-column figure
cairo_pdf(file = paste0(output.dir, "/junk.pdf"),
          width = 180/25.4, height = 150/25.4, onefile = TRUE)




# routine to run the process on the 
diff.tests <- list()
input_filename_list <- list.files(path = campaign_dir, pattern = campaign.filename.pattern)
for(filenumber in 1:length(input_filename_list)){
    # drop the extension (if any) from the filename to get the root used for naming new outputs
    # note: special handling of indexing at length == 1 - does not work for all dot filenames
    filename_root <- unlist(strsplit(basename(input_filename_list[filenumber]), split = "[.]"))
    filename_root <- paste(filename_root[1:length(filename_root)-1], collapse = ".")

    # create the directory for output if not already existing - parents too
    # remove spaces and underscores as per http://yihui.org/2018/03/space-pain
    compare_dir_name <- paste0(output.dir, "/", gsub("[ _]", "-", filename_root), "-compare-output")
    if (!dir.exists(compare_dir_name)){
        dir.create(compare_dir_name, recursive = TRUE)
    }

    # read the campaign file
    testing <- read.table(file.path(campaign_dir, input_filename_list[filenumber]),
                          sep = ",", quote = "\"", header = TRUE, stringsAsFactors = FALSE)

    # check for the expected columns and skip any file without
    if(length(unique(colnames(testing)[colnames(testing) %in%
               c("Type", "Filename", "Columns", "Dictionary", "Labels", "OrderTable")])) == 6){

    print(paste0("Processing: ", filename_root))
    testing$Type <- factor(testing$Type)
    testing$Columns <- lapply(strsplit(testing$Columns, ","), as.integer)

    # R converts blank strings into NA on input so reverse this before extracting the list
    testing$ExcludedSDBs[is.na(testing$ExcludedSDBs)] <- ""
    testing$ExcludedSDBs <- lapply(strsplit(testing$ExcludedSDBs, ","), trimws)

    testing_data <- read_and_merge_columns(testing, LiGA.data.dir = LiGA.data.dir,
            dict_dir = "dictionaries", order_table_dir = order_table_dir)
    if(NROW(testing_data$Contrasts) < 1) {
        print("No Test - Control contrasts available.")
    }

    #sum count table
    sum.table <- data.frame(Glycan = testing_data$table$Barcode,
                            SDB = testing_data$table$SDB,
                            testing_data$table[,3:(ncol(testing_data$table)-1)],
                            check.names = FALSE)
    # delete unkown SDBs
    ###sum.table <- sum.table[sum.table$Glycan != "",]

    # write the count table
    write.table(sum.table,
                file = paste0(compare_dir_name, "/", "count-table-",
                              sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)



# names of the internal standards
is.SDBs <- c("SDB267", "SDB271", "SDB272", "SDB273", "SDB274", "SDB275", "SDB276",  "SDB277", "SDB278", "SDB280")

    is.indices <- which(testing_data$table$SDB %in% is.SDBs)
    count.columns <- c(3:(ncol(testing_data$table)-1))
    ### # move all the IS counts to the first IS row and delete other rows
    is.summed <- colSums(testing_data$table[is.indices,count.columns])
    ### testing_data$table[is.indices[1],count.columns] <- is.summed
    ### testing_data$table <- testing_data$table[-is.indices[2:length(is.indices)],]

    # add a row to the count table and dictionary for the summed IS SDBs
    testing_data$table[nrow(testing_data$table) + 1,count.columns] <- is.summed 
    testing_data$table$SDB[nrow(testing_data$table)] <- "AzOH"
    testing_data$table$Barcode[nrow(testing_data$table)] <- "AzOH"
    testing_data$table$Order[nrow(testing_data$table)] <- 0
    testing_data$table <- testing_data$table[order(testing_data$table$Order),]
    testing_data$dict[nrow(testing_data$dict) + 1,1:2] <- c("AzOH", "AzOH")
    testing_data$dict$Order[nrow(testing_data$dict)] <-
                                      min(testing_data$dict$Order)-1
    is.SDBs <- c("AzOH")

#    # patch the contrasts to repeat tests against Naive
#    if("Naive" %in% colnames(testing_data$Column.types)){
#        non.naive.tests <- testing_data$Contrasts[testing_data$Contrasts$Reference != "Naive","Type"]
#        testing_data$Contrasts <- rbind(testing_data$Contrasts,
#            data.frame(Type = unique(non.naive.tests), Reference = "Naive"))
#    }

    # loop over three subsets of the data
    for(case in 1:1){
        if(case == 1){
            analysis_case <- "all"
            active_rows <- TRUE  #!is.na(testing_data$sequence)
        } 
        if(case == 2){
            analysis_case <- "knownSDBs"
            active_rows <- !is.na(testing_data$table$SDB)
        } 
        if(case == 3){
            analysis_case <- "knownGlycans"
            active_rows <- testing_data$table$Barcode != ""
        } 

        # use edgeR, apply TMM and calculate the mean & st. dev. for each condition

        # generate a mask to select data rows/columns  
        count.mask <- apply(testing_data$Column.types, 1, any)

        # generate a placeholder group - select the first if multiple groups
        a.group <- apply(testing_data$Column.types, 1,
                             function(x){head(names(which(x)), n=1)})

        dge <- DGEList(counts = testing_data$table[active_rows,count.mask],
                      group = factor(unlist(a.group[count.mask])))
        dge <- calcNormFactors(dge, method = "TMM")
       
        
        # normalize a "library" of just the internal standards
        # is.SDBs - names of the internal standards - set further up
        is.rows <- (testing_data$table$SDB %in% is.SDBs) & active_rows
        is.dge <- DGEList(counts = testing_data$table[is.rows,count.mask],
                          group = factor(unlist(a.group[count.mask])))
        is.dge$samples$lib.size <- dge$samples$lib.size
        is.dge <- calcNormFactors(is.dge, method = "RLE")  # I think RLE is most appropriate uses all rows
        
        # copy the normalization factors to the whole library and calculate CPM
        dge$samples$norm.factors <- is.dge$samples$norm.factors
        dge.normalized <- cpm(dge, normalized.lib.sizes = TRUE)

### # only use row #10 to normalize TMM
### a <- testing_data$table[active_rows, count.mask][10,]
### b <- apply(testing_data$table[active_rows, count.mask], 2, sum)
### c <- rbind(a, a, a, a, b - 4*a)
### dge_r <- DGEList(counts = c, group = factor(unlist(a.group[count.mask])))
### dge_r <- calcNormFactors(dge_r, method = "TMM")
### dge$samples$norm.factors <- dge_r$samples$norm.factors
        
        # extract the subset of the test and control data
        testing_data$Contrasts$Type      <- as.character(testing_data$Contrasts$Type)
        testing_data$Contrasts$Reference <- as.character(testing_data$Contrasts$Reference)
        group.names <- data.frame(Group = unique(c(testing_data$Contrasts$Type,
                                                   testing_data$Contrasts$Reference)), stringsAsFactors = FALSE)
        group.names$Name <- make.names(unique(group.names$Group), unique = TRUE, allow_ = FALSE)
        contrast.table <- left_join(left_join(testing_data$Contrasts, group.names, by = c("Type" = "Group")),
                                                                      group.names, by = c("Reference" = "Group"))
        contrast.table$String <- paste0(contrast.table$Name.x, "-", contrast.table$Name.y)

        #TODO: Would it be better to estimate dispersion from ALL the data?
        for(pair in 1:nrow(contrast.table)){
            # extract masks selecting which datasets are used in the current contrast        
            type.data       <- testing_data$Column.types[,contrast.table$Type[pair]]
            compare.to.data <- testing_data$Column.types[,contrast.table$Reference[pair]]

            # need omit any datasets on both sides of the comparison    
            both.data <- type.data & compare.to.data
            if(any(both.data)){
                stop("Dataset(s) used on both sides of comparison:", colnames(testing_data$table)[both.data])
            }             

            # combine the indexing masks and extract the group names
            subset_indices <- (type.data | compare.to.data) & !both.data
            subset_group <- ifelse(testing_data$Column.types[subset_indices,contrast.table$Type[pair]],
                    contrast.table$Type[pair], contrast.table$Reference[pair])

            # edgeR object for the subset, take norm.factors from whole TMM
            dge_subset <- DGEList(counts = testing_data$table[active_rows,subset_indices],
                                 norm.factors = dge$samples$norm.factors[subset_indices[count.mask]],
                                 group = subset_group)

            # identify the low-count data (rows) in the subset
            Valid <- filterByExpr(dge_subset)
            #dge_subset <- dge_subset[Valid, , keep.lib.sizes = FALSE]

            # the following is from the tutorial / Jessica's script, needs to be further modified
            # build the model matrix and fix-up the output names (meet R's variable name rules)
            names <- left_join(data.frame(Group = as.character(dge_subset$samples$group),
                                 stringsAsFactors = FALSE), group.names, by = "Group")$Name
            design <- model.matrix(~ 0 + names, data = dge_subset$samples)
            colnames(design) <- gsub("^names", "", colnames(design))

            # estimate dispersion (by one method or another)
            tmp <- try(estimateDisp(dge_subset, design, robust = TRUE), TRUE)
            if(inherits(tmp, "try-error")){
                dge_subset <- estimateCommonDisp(dge_subset)
                dge_subset <- estimateTagwiseDisp(dge_subset)
            } else {
                # print("Using estimateDisp")
                dge_subset <- estimateDisp(dge_subset, design, robust = TRUE)
                ###plotMDS(dge_subset, top = dim(dge_subset$counts)[1], col = colgroup)
            }

            # now use GLM to model the glycan dispersion and test for differences
            fit <- glmQLFit(dge_subset, design, robust = TRUE)
            contrasts <- makeContrasts(contrasts = contrast.table$String[pair], levels = design)
            qlfi <- glmQLFTest(fit, contrast = contrasts)

            # build a data.frame of output information...
            difference.frame <- data.frame(
                    rank = 1:sum(active_rows),
                    type.rank = 1 + length(diff.tests),
                    testing_data$table[active_rows, c("SDB", "Barcode", "Order")],
                    contrast.table[pair,1:2])

            # ...include the TMM normalized mean and st. deviation...
            # to calculate relative error for the error bars
            difference.frame[["Test_CPM"]]    <- apply(dge.normalized[,type.data[count.mask]], 1, mean)
            difference.frame[["Test_SD"]]     <- apply(dge.normalized[,type.data[count.mask]], 1, sd)
            difference.frame[["Control_CPM"]] <- apply(dge.normalized[,compare.to.data[count.mask]], 1, mean)
            difference.frame[["Control_SD"]]  <- apply(dge.normalized[,compare.to.data[count.mask]], 1, sd)

            # join the quasi-liklihood F-test results back the appropriate IDs
            # NAs are inserted for low-count cases excluded from the testing
            
            difference.frame <- left_join(difference.frame,
                data.frame(testing_data$table[active_rows, c("SDB", "Barcode", "Order")],
                    Test.CPM = exp(fit$coefficients[,contrast.table$Name.x[pair]])*1e6,
                    Control.CPM = exp(fit$coefficients[,contrast.table$Name.y[pair]])*1e6,
                    qlfi$table, QValue = p.adjust(qlfi$table$PValue, "BH"), Valid),
                by = c("SDB", "Barcode", "Order"))
            #difference.frame$PValue[!Valid] <- NA

            bar.plot <- make.barplot(difference.frame, TRUE, filename_root,
                                     font.family = font.family)

            diff.tests[[length(diff.tests)+1]] <- list(frame = difference.frame,
                      plot = bar.plot,
                      data1.set = type.data,
                    data2.set = compare.to.data)


#            render("../report-revised.Rmd", ###output_format = "pdf_document",
#                   params = list(testing_data = testing_data,
#                                 difference.frame = difference.frame,
#                                 bar.plot = bar.plot[[1]]))


            compare.table <- difference.frame[order(difference.frame$PValue),
                c("Barcode", "SDB", "Test.CPM", "Test_SD", "Control.CPM",  
                  "Control_SD", "logFC", "QValue")] 
            colnames(compare.table)[3:8] <- c("Test (CPM)", "Test (SD)",
                    "Reference (CPM)", "Ref (SD)", "FC", "q-value")
            compare.table$FC <- 2**compare.table$FC
            

###            compare.table <- cbind(compare.table,
###                              sum.table[match(compare.table$SDB, testing_data$table$SDB),
###                                        3:NCOL(testing_data$table)])

            row.idx <- match(testing_data$table$SDB, compare.table$SDB)
            row.flags <- !is.na(row.idx)
###            ratmir.table <- cbind(testing_data$table[row.flags,-ncol(testing_data$table)],
###                      compare.table[row.idx[row.flags],])


            output.order <- order(testing_data$table$Order)
            ratmir.table <- cbind(compare.table[row.idx[row.flags],], 
                   testing_data$table[output.order, type.data],
                   testing_data$table[output.order, compare.to.data]) 
            write.table(ratmir.table,
                file = paste0(compare_dir_name, "/", "data-table-",
                      sub(".[^.]*$", "", input_filename_list[filenumber]), "-",
                      contrast.table$String[pair], ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
         }

#        # sanitize the column names and values (remove tabs) and write the table
#        output_table <- output_df
#        colnames(output_table) <- gsub("\t", " ", colnames(output_table))
#        output_table$Barcode <- gsub("\t", " ", output_table$Barcode)
#        output_table$SDB     <- gsub("\t", " ", output_table$SDB)
#        data_table_filename2 <- file.path(".", paste0(filename_root, "+", analysis_case, ".tsv"))
#        write.table(output_table,
#                    file = data_table_filename2, sep = "\t", quote = FALSE, row.names = FALSE)




    }
    }

}






# arrange the difference tests into one large table for plotting
long.data <- do.call(rbind, lapply(diff.tests, function(x){x$frame}))
#long.data$PValue <- ifelse(long.data$Valid, long.data$PValue, NA)
long.data$QValue <- p.adjust(long.data$PValue, method = "BH")

long.data$Glycan <- sub("-\\[[0-9<]*]$", "",  long.data$Barcode)
long.data$Density <- suppressWarnings(as.integer(
                          sub(".*-\\[<?0*([1-9][0-9]*)\\]$", "\\1", long.data$Barcode)))

permutation <- order(long.data$Order,			# primary sort by provided Order
                     long.data$Barcode == "",           # empty labels at the end
                     long.data$Glycan,		        # group each glcan's entries (alphabetic)
                     long.data$Density,
                     #long.data$Barcode,	                # finally just use the raw string
                     nchar(long.data$SDB),		# resolve ties by the SDB, but put
                     long.data$SDB, na.last = TRUE)	# shorter ones first (numeric order)

long.data <- long.data[permutation,]
x.position.frame <- data.frame(unique(long.data[,c("Order", "Barcode", "SDB")]), stringsAsFactors = FALSE)
x.position.frame$xpos = 1:nrow(x.position.frame)
long.data <- left_join(long.data, x.position.frame, by = c("SDB", "Barcode", "Order"))

long.data <- left_join(long.data, data.frame(type.rank = c(1:7), ypos=c(7, 4, 3, 2, 5, 6, 1)), by = "type.rank")


        # sanitize the column names and text values (remove tabs) and write the table
        output_table <- long.data[,c("Barcode", "SDB", "Type", "Reference",
            "Test.CPM", "Test_SD", "Control.CPM", "Control_SD",
            "logFC", "F", "PValue", "QValue")]
        colnames(output_table) <- gsub("\t", " ", colnames(output_table))
        output_table$Barcode   <- gsub("\t", " ", output_table$Barcode)
        output_table$SDB       <- gsub("\t", " ", output_table$SDB)
        output_table$Type      <- gsub("\t", " ", output_table$Type)
        output_table$Reference <- gsub("\t", " ", output_table$Reference)
        data_table_filename2 <- file.path(compare_dir_name,
                paste0(filename_root, "+", analysis_case, ".tsv"))
        ###write.table(output_table,
        ###            file = data_table_filename2, sep = "\t", quote = FALSE, row.names = FALSE)

long.data$Label <- paste0(sub("'", "\u2032", long.data$Barcode), " ", long.data$SDB)
long.data$Label <- reorder(factor(long.data$Label), long.data$xpos)

significant.rows.1 <- long.data$Valid & long.data$QValue <= 0.05 & log2(2) <= long.data$logFC & long.data$logFC < log2(3.8)
significant.rows.2 <- long.data$Valid & long.data$QValue <= 0.05 & log2(3.8) <= long.data$logFC

#sum count table
sum.table <- data.frame(Glycan = testing_data$table$Barcode,
                       SDB = testing_data$table$SDB,
                       testing_data$table[,3:(ncol(testing_data$table)-1)],
                       check.names = FALSE)

# delete unkown SDBs
###sum.table <- sum.table[sum.table$Glycan != "",]

# write the count table
write.table(sum.table,
        file = paste0(compare_dir_name, "/", "count-table-2nd-", sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE)

#comparison table
compare.table <- long.data[order(long.data$PValue), c("Barcode", "SDB", "Type", "Reference",
             "Test.CPM", "Test_SD", "Control.CPM", "Control_SD", "logFC", "F", "PValue", "QValue")]
colnames(compare.table)[3:12] <- c("Test", "Reference", "Test (CPM)", "Test (SD)",
                                   "Reference (CPM)", "Ref (SD)", "log2(FC)", "F", "p-value", "q-value")

compare.table <- cbind(compare.table, sum.table[match(compare.table$SDB, sum.table$SDB),3:NCOL(sum.table)])
write.table(compare.table,
        file = paste0(compare_dir_name, "/", "comparison-table-", sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE)


row.idx <- match(testing_data$table$SDB, compare.table$SDB)
row.flags <- !is.na(row.idx)
ratmir.table <- cbind(testing_data$table[row.flags,-ncol(testing_data$table)],
                      compare.table[row.idx[row.flags], c("Test (CPM)", "Test (SD)", "Reference (CPM)", "Ref (SD)", "log2(FC)", "q-value")])
colnames(ratmir.table)[ncol(ratmir.table)-1] <- "FC"
ratmir.table$FC <- 2**ratmir.table$FC
write.table(ratmir.table,
        file = paste0(compare_dir_name, "/", "data-table-", sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE)



###
###
###
###
        plots <- list()
        for(i in 1:length(rownames(dge$samples))){
             row.name <- rownames(dge$samples)[i]
             group.name <- dge$samples$group[i]
             raw.column <- compare.table[row.name][,1]   # second index strips colname
             if(group.name == "Test"){
                 ref.column <- compare.table$"Test (CPM)"
             } else {
                 ref.column <- compare.table$"Reference (CPM)"
             }
             raw.mask <- raw.column > 0    # remove zero counts for log-scale
             subplot <- ggplot(data.frame(x = ref.column[raw.mask], y = raw.column[raw.mask])) +
                            geom_point(aes(x = x, y = y), size = 1) + scale_x_log10() + scale_y_log10() +
                            labs(x = row.name, y = "") +
                            theme(text = element_text(family = font.family, size = 7))
             #print(subplot)
             plots[[i]] <- ggplotGrob(subplot)
        }
        ###do.call(grid.arrange, plots)
#        for(i in 1:ceiling(length(plots)/12)){
#            do.call(grid.arrange, c(plots[(12*i-11):min(12*i,length(plots))],
#                                    ncol = 4, nrow = 3))
#        }








long.data$Type <- factor(long.data$Type, levels = unique(long.data$Type))

##### special handling for the data in question - EC 2022-03-23
long.data$yPos <- -2/3 * c(0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3) - as.numeric(long.data$Type)
###########


### replace to CD22 cell data with the logFC from:
### Feb/outputs/Dec-CD22-cell-SC-4-compare-output/comparison-table-Dec-CD22-cell\ SC-4.tsv

custom.logFC <- c(1.82203288557042, -0.682016620146907, 8.10107800362924, 10.3322836817092, 0, 8.01843130294255, 1.84349264695734, 0.873390317917666, -5.912631508891, 5.06983568159416, 11.3061974186283, 5.35253646635049, -3.64176235173184, 3.21639221883177, 10.3013277928155, 5.24900344454714, 10.5112154611616, 5.74228453208128, 0, -1.42086406971501, 2.55101436541334, -1.71575307611652, 0, -2.57822450039994, 0.255411422988217, -2.87824783622648, 2.62596375253971, -1.50247148262152, -4.34402424960869, -0.6622502955705, -0.296666064936306, -7.64378672072682, -7.26304579096009, -7.15725598527302, -2.14837252077867, -6.96048984955859, -1.04325829020026, 0.874805469411235, -3.68656837097469)

custom.SDB <- c("SDB103", "SDB191", "SDB28", "SDB120", "SDB3", "SDB186", "SDB36", "SDB206", "SDB163", "SDB170", "SDB43", "SDB171", "SDB108", "SDB219", "SDB110", "SDB158", "SDB156", "SDB147", "SDB41", "SDB101", "SDB119", "SDB37", "SDB93", "SDB13", "SDB221", "SDB90", "SDB102", "SDB155", "SDB94", "SDB29", "SDB222", "SDB224", "SDB104", "SDB136", "SDB95", "SDB20", "SDB225", "SDB17", "SDB267")

###custom.SDB[custom.SDB == "SDB41"] <- "SDB139"  # both glycan# 10-[50] barcodes
indices <- match(long.data[long.data$Type == "CD22 cell",]$SDB, custom.SDB)
long.data[long.data$Type == "CD22 cell",]$logFC <- custom.logFC[indices]

# restrict the heatmsap the expected labels
long.data <- long.data[long.data$Barcode != "" & !grepl("^AzOH", long.data$Barcode),]


# recalculate the FDR on the actual dataset
long.data$QValue <- p.adjust(long.data$PValue, method = "BH")

vert.range = c(min(long.data$yPos)-0.5, max(long.data$yPos)+0.5)
significant.rows <- long.data$Valid & long.data$QValue <= 0.05 & long.data$logFC > 1


heatmap <- ggplot(long.data,
                  aes(x = Label, y = yPos)) +
  theme_light() +
  theme(text = element_text(family = font.family, size = 7, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  coord_fixed(ratio = 0.653) +
  

  geom_tile(aes(fill = pmax(0.1, 2 ** logFC)), colour = "gray70", size = 0.2343) +
  geom_vline(xintercept = 0.5 + c(5, 11, 17, 23, 29), colour = "black", size = 0.4686) +
  scale_y_continuous(breaks = unique(long.data$yPos),
                     labels = levels(long.data$Type),
                     sec.axis = dup_axis(breaks = c(-20/3, -2, -71/6, -16.5),
                     labels = c("Mannose\nBinders", "Neu5Ac\nBinders",
                                "Gal\nBinders", "GlcNAc\nBinders"))) +

  geom_text(label = "*", hjust = 0.45, vjust = 0.75, colour = "white", 
            data = long.data[significant.rows,]) +
  labs(x = "Glycan", y =  "Binding Target", aes(hjust = 0)) +
  theme(axis.title = element_text(size = 7.0),
        axis.text = element_text(colour= "black", size = 7.0),
        axis.text.x  = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 7.0),
        axis.text.y  = element_text(hjust = 0, size = 7.0), 
        axis.title.y.right = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(), # element_text(size = 10),
        legend.text = element_text(size = 7.0),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(2.5, "cm"))

cb.max <- 2550
colour_steps <- colour_gradient()
#colour_steps <- colour_steps[54:256]  reduce the range of the colour scale
colour_intervals <- 0:(length(colour_steps) - 1) / (length(colour_steps) - 1)


colour.intervals.B <- (log(100 * colour_intervals^2) - log(0.1)) / (log(cb.max) - log(0.1))
extra.points <- (log(c(500, cb.max)) - log(0.1)) / (log(cb.max) - log(0.1))
colour.intervals.B <- c(colour.intervals.B, extra.points)
colour.steps.B <- c(colour_steps, "#ff0000", "#ffb0b0")


# height and width are in inches - width choosen for a two-column figure
cairo_pdf(file = paste0(compare_dir_name, "/", heatmap.filename),
          width = 180/25.4, height = 150/25.4, onefile = TRUE)


print(heatmap + scale_fill_gradientn(colours = colour.steps.B, values = colour.intervals.B,
                       limits = c(0.1, cb.max), trans = "log",
                       breaks = c(0.1, 1, 10, 100, 1000, 2500),
                       labels = c("0.1", "1", "10", "100", "1000", "2500"),
                       na.value = "white",  guide = guide_colorbar(nbin = 1000))
)

# end the PDF file generation
dev.off()


# use knitr to produce the PDF report
render("report.Rmd", output_dir = compare_dir_name, intermediates_dir = compare_dir_name,
	###output_format = "pdf_document",
       params = list(testing_data = testing_data,
                     difference.tests = diff.tests,
                     main.table = long.data))
#                     bar.plot = bar.plot[[1]]))

dev.off()
