#' Pull and organise chromosome ideogram data from HGSC
#'
#' @param chromosome Chromosome to display by number. Typically autosomes are listed first in numerical order, then non-homologous chromosomes.
#' @param genome Abbreviated genome name for data to be pulled
#' @param Y_offset y-coordinate for the centre of the ideogram
#' @param height Height of the ideogram
#'
#' @export
#'

geom_ideogram <- function(chromosome, genome="hg38", Y_offset=0, height=1){
  chromosome_data <- generate_drawing_data(genome=genome)
  chrs=stringr::str_sort(unique(chromosome_data$arms$Chromosome), numeric = TRUE)
  list(
    ggplot2::layer(
      geom = "rect",
      mapping = ggplot2::aes(xmin=Start,xmax=End,ymin=Y_offset-(height/2), ymax=Y_offset+(height/2)),
      data = chromosome_data$arms[chromosome_data$arms$Chromosome==chrs[chromosome],],
      params = list(fill=chromosome_data$arms[chromosome_data$arms$Chromosome==chrs[chromosome],]$Color),
      stat = "identity",
      position = "identity",
      show.legend = FALSE,
      inherit.aes = FALSE),
    ggplot2::layer(
      geom = "polygon",
      mapping = ggplot2::aes(X,Y_offset - Y * (height / 2)),
      data=chromosome_data$centromeres[chromosome_data$centromeres$Chromosome==chrs[chromosome],],
      params=list(fill=chromosome_data$centromeres[chromosome_data$centromeres$Chromosome==chrs[chromosome],]$Color),
      stat = "identity",
      position = "identity",
      show.legend = FALSE,
      inherit.aes = FALSE),
    ggplot2::layer(
      geom = "path",
      mapping = ggplot2::aes(X,Y_offset - Y * (height / 2)),
      data = chromosome_data$outlines[chromosome_data$outlines$Chromosome==chrs[chromosome],],
      params = list(color="black"),
      stat = "identity",
      position = "identity",
      show.legend = FALSE,
      inherit.aes = FALSE)
  )
}

get_ideo_data <- function(genome="hg38", strip=TRUE){
  down_text <- readLines(gzcon(url(paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",genome,"/database/cytoBand.txt.gz"))))
  cyto <- utils::read.csv(textConnection(down_text), sep="\t", header=FALSE)
  colnames(cyto) <- c("Chromosome","Start","End","Name","Type")
  if (strip){
    cyto <- cyto[! grepl("_",cyto$Chromosome),]
  }
  cyto$Color <- "white"
  cyto$Color[cyto$Type=="gpos25"] <- "grey25"
  cyto$Color[cyto$Type=="gpos50"] <- "grey50"
  cyto$Color[cyto$Type=="gpos75"] <- "grey75"
  cyto$Color[cyto$Type=="gpos100"] <- "black"
  cyto$Color[cyto$Type=="gvar"] <- "powderblue"
  cyto$Color[cyto$Type=="acen"] <- "indianred1"

  for (chr in unique(cyto$Chromosome)){
    if (dim(cyto[cyto$Chromosome==chr & cyto$Type=="acen",])[1]==0){
      cyto <- cyto[cyto$Chromosome != chr,]
    }
  }

  return(cyto)
}


generate_drawing_data <- function(genome="hg38"){
  centromeres <- data.frame("Chromosome"=NA,"X"=NA,"Y"=NA,"Type"=NA,"Color"=NA)
  cyto <- get_ideo_data(genome=genome, strip=TRUE)
  chrs <- unique(cyto$Chromosome)
  for (chr in chrs){
    centromeres <-rbind(centromeres,
                        data.frame("Chromosome"=chr,"X"=cyto[cyto$Chromosome==chr & cyto$Type=="acen","Start"], "Y"=c(1,0.6),"Type"="acen","Color"="indianred1"),
                        data.frame("Chromosome"=chr,"X"=max(cyto[cyto$Chromosome==chr & cyto$Type=="acen","End"]), "Y"=c(1,-1),"Type"="acen","Color"="indianred1"),
                        data.frame("Chromosome"=chr,"X"=sort(cyto[cyto$Chromosome==chr & cyto$Type=="acen","Start"],decreasing = TRUE), "Y"=c(-0.6,-1),"Type"="acen","Color"="indianred1"))
  }
  centromeres <- stats::na.omit(centromeres)
  chromosome_outlines <- data.frame("Chromosome"=NA,"X"=NA,"Y"=NA)
  for (chr in chrs){
    chromosome_outlines <- rbind(chromosome_outlines,
                                 data.frame("Chromosome"=chr, "X"=c(0, centromeres[centromeres$Chromosome==chr,"X"][1:3], max(cyto[cyto$Chromosome==chr,"End"]),max(cyto[cyto$Chromosome==chr,"End"]),centromeres[centromeres$Chromosome==chr,"X"][4:6],0,0),"Y"=c(1,1,0.6,1,1,-1,-1,-0.6,-1,-1,1)))
  }
  chromosome_outlines <- stats::na.omit(chromosome_outlines)
  return(list(arms=cyto[cyto$Type != "acen",], centromeres=centromeres, outlines=chromosome_outlines))
}






