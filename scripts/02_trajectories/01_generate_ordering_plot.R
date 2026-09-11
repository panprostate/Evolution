#tumour_set = "SherlockLung_test_bundle55SD_DN1"
#input_dir = "Sherlock/New_Histology/2025/Data/Test_Bundle/"
#driver_file = "Sherlock/sherlock_driver_symbols.txt"
#cyto_file = "Sherlock/cytoBand_hg38.txt"
#assembly = "hg38"
#PL_method = "PlackettLuce"
#wgd_file = "Sherlock/SherlockLung_WGDStatus.txt"
#including_wgd = TRUE
#including_mrca = FALSE

tumour_set = ""

input_dir = paste0("outputs/02_trajectories/", tumour_set, "/")
driver_file = "data/raw/Somatic_drivers/DRIVER_ANNOTATIONS__2023-09-25_release4__SNV+indel__filteredVariants__driverGenes.txt"
cyto_file = "data/ext/cytoBand_hg19.txt"
assembly = "hg19"
PL_method = "PlackettLuce"
wgd_file = "data/processed/wgd/PPCG_WGD_status_PCAWG_method_20260212.txt"
including_wgd = FALSE
including_mrca = TRUE
use_remote_db = FALSE

output_dir = input_dir
figure_dir = "figures/02_trajectories/"
#output_dir = "Sherlock/New_Histology/2025/Figures/Fig1/"
#output_dir = "Sherlock/New_Histology/2025/Figures/Obfuscated/"
horizontal = FALSE

cat(paste0("Running with the following variable values:

tumour_set = ", tumour_set, "
input_dir = ", input_dir, "
driver_file = ", driver_file, "
cyto_file = ", cyto_file, "
assembly = ", assembly, "
PL_method = ", PL_method, "
wgd_file = ", wgd_file, "
including_wgd = ", including_wgd, "
including_mrca = ", including_mrca, "

"))

### INITIALISATION ###

library(biomaRt)
library(GenomicRanges)
library(magrittr)
library(hrbrthemes)
# library(extrafont)
library(ggplot2)
library(gtools)

geom_stripes <- function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         ...,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomStripes,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(...)
  )
}

GeomStripes <- ggplot2::ggproto("GeomStripes", ggplot2::Geom,
                                required_aes = c("y"),
                                
                                default_aes = ggplot2::aes(
                                  ymin = -Inf, ymax = Inf,
                                  #odd = "#22222222", even = "#00000000",
                                  odd = "#99999911", even = "#FFFFFF00",
                                  # Change 'size' below from 0 to NA.
                                  # When not NA then when *printing in pdf device* borders are there despite
                                  # requested 0th size. Seems to be some ggplot2 bug caused by grid overriding
                                  # an lwd parameter somewhere, unless the size is set to NA. Found solution here
                                  # https://stackoverflow.com/questions/43417514/getting-rid-of-border-in-pdf-output-for-geom-label-for-ggplot2-in-r
                                  alpha = NA, colour = "black", linetype = "solid", size = NA
                                ),
                                
                                # draw_key = ggplot2::draw_key_blank,
                                draw_key = ggplot2::draw_key_rect,
                                
                                draw_panel = function(data, panel_params, coord) {
                                  ggplot2::GeomRect$draw_panel(
                                    data %>%
                                      dplyr::mutate(
                                        x = round(.data$x),
                                        xmin = .data$x - 0.5,
                                        xmax = .data$x + 0.5
                                      ) %>%
                                      dplyr::select(
                                        .data$ymin, .data$ymax,
                                        .data$xmin, .data$xmax,
                                        .data$odd, .data$even,
                                        .data$alpha, .data$colour, .data$linetype, .data$size
                                      ) %>%
                                      unique() %>%
                                      dplyr::arrange(.data$xmin) %>%
                                      dplyr::mutate(
                                        .n = dplyr::row_number(),
                                        fill = dplyr::if_else(
                                          .data$.n %% 2L == 1L,
                                          true = .data$odd,
                                          false = .data$even
                                        )
                                      ) %>%
                                      dplyr::select(-.data$.n, -.data$odd, -.data$even),
                                    panel_params,
                                    coord
                                  )
                                }
)


cyto = read.table(cyto_file, stringsAsFactors=FALSE, header=FALSE)
cyto_gr = GRanges(seqnames=gsub("chr", "", cyto[,1]), ranges=IRanges(start=cyto[,2], end=cyto[,3]), cytoband=cyto[,4])

### GENERATE PLOTS ###

G = 1
trajectories = 1:3

for (trajectory in trajectories){
  message(paste0("Plotting time of events in trajectory ", trajectory, "..."))
  mergedseg=read.table(paste0(input_dir, "/", "PPCG_Feb2026_DN", trajectory, "_mergedseg_with_clonality.txt"),header=T,stringsAsFactors = F)
  matrix_coresp=read.table(paste0(input_dir, "/", "PPCG_Feb2026_DN", trajectory, "_matrix_coresp_G1.txt"),header=T,stringsAsFactors = F)
  matrix_out=read.table(paste0(input_dir, "/", "PPCG_Feb2026_DN", trajectory, "_matrix_out_G1.txt"),header=T,stringsAsFactors = F)
  mergedseg$ID = gsub("chr23", "chrX", mergedseg$ID)

  cnai = which(mergedseg$CNA %in% c("dGain","dLOH","dHD"))
  mergedseg$ID[cnai] = sapply(cnai, function(i) { paste0(gsub("d", "", mergedseg$CNA[i]), "_", mergedseg$ID[i]) })

  if (wgd_file != "none") {
    wgd_status = read.table(wgd_file, stringsAsFactors=FALSE)
    all_wgd_samples = wgd_status[wgd_status[,2] == "WGD",1]
    wgd_samples = list(intersect(unique(mergedseg$Tumour_Name), all_wgd_samples))
  } else {
    wgd_samples = list(unique(mergedseg[mergedseg$tumour_ploidy > 2,'Tumour_Name']))
  }

  if (length(wgd_samples[[1]]) == 0) {
    including_wgd = FALSE
    print("No WGD samples found - will not attempt to include WGD")
  }

  event.id=unique(mergedseg[,c("ID","noevent","CNA")])
  if (including_mrca) {
    event.id=rbind(data.frame(ID="MRCA",noevent=(length(unique(mergedseg$ID))+1),CNA="cMRCA",stringsAsFactors = F),event.id)
  }
  if (including_wgd) {
    event.id=rbind(data.frame(ID="Whole_Genome_Duplication",noevent=(length(unique(mergedseg$ID))+1),CNA="cWGD",stringsAsFactors = F),event.id)
  }
  event.id$cna=base::substring(event.id$CNA,first = 2)
  event.id=unique(event.id[,c("ID","noevent","cna")])

  #event.id$ID = gsub("chr23", "chrX", event.id$ID)

  names(matrix_coresp)=c("noevent","event")
  if (PL_method == "PLMIX") {
    matrix_coresp$event=paste("p",matrix_coresp$event,sep = "_")
  }
  coresp=merge(event.id,matrix_coresp,"noevent")

  values <- reshape2::melt(matrix_out,id.vars=NULL)
  values$Var1=NULL
  names(values)=c("event","value")
  values$event = gsub("X","",values$event)

  out <- merge(values,coresp,"event")

  PLOTDF=data.frame()
  unoevent=unique(out$noevent)
  for (i in 1:length(unique(out$noevent))){
    event=out[out$noevent==unoevent[i],]
    event=event[order(event$value),]
    event=event[26:975,]
    PLOTDF=rbind(PLOTDF,event)
  }

  all_event_i = unique(PLOTDF$event)
  event_medians = sapply(all_event_i, function(event_i) {
    median(PLOTDF[PLOTDF$event == event_i,'value'])
  })
  PLOTDF$ID=factor(PLOTDF$ID,levels=unique(PLOTDF$ID)[order(event_medians, decreasing=T)])

  unique_IDs = unique(mergedseg$ID)
  samples_with_each_event = lapply(unique_IDs, function(ID) {
    unique(mergedseg[mergedseg$ID == ID,'Tumour_Name'])
  })
  if (including_wgd) {
    samples_with_each_event = append(samples_with_each_event, wgd_samples)
    names(samples_with_each_event) = c(unique_IDs, "WGD")
  } else {
    names(samples_with_each_event) = unique_IDs
  }
  unique_samples = unique(mergedseg$Tumour_Name)
  n_samples = length(unique_samples)
  if (including_mrca) {
    new_names = c(names(samples_with_each_event), "MRCA")
    samples_with_each_event = append(samples_with_each_event, list(unique_samples))
    names(samples_with_each_event) = new_names
  }
  event_frequency = sapply(samples_with_each_event, function(x) { round((length(x) / n_samples) * 100, 1) })
  names(event_frequency) = gsub("WGD", "Whole_Genome_Duplication", names(event_frequency))
  event_frequency = event_frequency[levels(PLOTDF$ID)]

  if (driver_file != "none") {
    all_drivers = unique(read.table(driver_file, header=T, sep='\t', stringsAsFactors=FALSE)$symbol)
  }

  if (use_remote_db) {
  if (driver_file != "none") {
    thisMart <- biomaRt::useEnsembl(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl", 
                        GRCh = 37)
    driver_df = biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                              mart = thisMart,
                              filters = "hgnc_symbol",
                              values = all_drivers)
  } else {
    driver_df = c()
  }
    driver_gr = GRanges(seqnames=driver_df[,1], ranges=IRanges(start=driver_df[,2], end=driver_df[,3]), symbol=driver_df[,4])
    save(driver_gr, file=paste0("data/ext/driver_gene_gr_", assembly, ".Rdata"))
  } else {
    if (assembly == "hg38") {
      load("data/ext/driver_gene_gr_hg38.Rdata")
    } else {
      load("data/ext/driver_gene_gr_hg19.Rdata")
    }
    # driver_gr = all_gene_gr[all_gene_gr$symbol %in% all_drivers]
  }

  levels(PLOTDF$ID) = sapply(levels(PLOTDF$ID), function(id) {
    text_to_add = ""
    #if (substr(id, 1, 3) == "chr") {
    if (length(grep("_chr", id)) > 0) {
      locus_split = strsplit(id, "_")[[1]]
      #locus_gr = GRanges(seqnames=gsub(".*chr", "", locus_split[1]), ranges=IRanges(start=as.numeric(locus_split[2]), end=as.numeric(locus_split[3])))
      locus_gr = GRanges(seqnames=gsub(".*chr", "", locus_split[2]), ranges=IRanges(start=as.numeric(locus_split[3]), end=as.numeric(locus_split[4])))
      ol = suppressWarnings(findOverlaps(locus_gr, driver_gr))
      if (length(ol) > 0) {
        text_to_add = paste0(" (", paste(driver_gr$symbol[subjectHits(ol)], collapse=", "), ")")
      }
    }
    paste0(id, text_to_add)
  })
  levels(PLOTDF$ID) = gsub("dMut_", "", levels(PLOTDF$ID))
  levels(PLOTDF$ID) = gsub("Whole_Genome_Duplication", "WGD", levels(PLOTDF$ID))
  names(event_frequency) = levels(PLOTDF$ID)

  if (including_wgd) {
    WGD_linepos = -log10(sort(event_medians, decreasing=T)[which(levels(PLOTDF$ID) == "WGD")])
  }
  levels(PLOTDF$ID) = paste0(levels(PLOTDF$ID), " [", format(event_frequency, nsmall=1), "%]")

  cna_id_is = grep("_chr", levels(PLOTDF$ID))
  ids = sapply(strsplit(levels(PLOTDF$ID)[cna_id_is], " "), function(x) { paste(strsplit(x[1], "_")[[1]][2:4], collapse="_") })
  seg_mat = Reduce(rbind, strsplit(ids, "_"))
  seg_gr = GRanges(seqnames=gsub("chr", "", seg_mat[,1]), ranges=IRanges(start=as.numeric(seg_mat[,2]), end=as.numeric(seg_mat[,3])), id=ids)
  chr_cbs = lapply(split(cyto, cyto[,1]), function(x) {
    list(`p`=x[grep("p", x[,4]),4], `q`=x[grep("q", x[,4]),4])
  })
  names(chr_cbs) = gsub("chr", "", names(chr_cbs))
  ov = findOverlaps(seg_gr, cyto_gr)
  cyto_ids = sapply(unique(queryHits(ov)), function(qhi) {
    cbs = mcols(cyto_gr)[subjectHits(ov)[queryHits(ov) == qhi],'cytoband']
    chr = as.vector(seqnames(cyto_gr)[subjectHits(ov)[queryHits(ov) == qhi]])[1]
    if (length(cbs) > 1) {
      arm = substr(cbs[1], 1, 1)
      if (cbs[1] == chr_cbs[[chr]][[arm]][1] && cbs[length(cbs)] == chr_cbs[[chr]][[arm]][length(chr_cbs[[chr]][[arm]])]) {
        paste0(chr, arm)
      } else {
        #paste0(chr, cbs[1], "-", chr, cbs[length(cbs)])
        last_cb = cbs[length(cbs)]
        last_cb = substr(last_cb, 2, nchar(last_cb))
        paste0(chr, cbs[1], "-", last_cb)
      }
    } else {
      paste0(chr, cbs[1])
    }
  })
  levels(PLOTDF$ID)[cna_id_is] = sapply(1:length(cna_id_is), function(cna_i) {
    cna_id_i = cna_id_is[cna_i]
    splitty = strsplit(levels(PLOTDF$ID)[cna_id_i], " ")[[1]]
    part1 = strsplit(splitty[1], "_")[[1]][1]
    part2 = cyto_ids[cna_i]
    part3 = paste(splitty[2:length(splitty)], collapse=" ")
    paste0(part1, ": ", part2, " ", part3)
  })

  library(tidyverse)
  tdata = as_tibble(PLOTDF)

  #min_freq <- 3
  min_freq <- 7.5
  #min_freq <- 0
  tdata <- tdata %>% 
    arrange(ID) %>% 
    mutate(freq=str_remove(str_remove(ID,"^.*\\["),"%\\].*")) %>% 
    mutate(freq=as.numeric(freq)) %>% 
    filter(freq>min_freq)
    #filter(freq>=min_freq)

  saveRDS(tdata, file=paste0(output_dir, "/", trajectory, "_PLPlot_data.rds"))
  # filtering by minimum frequency may get rid of WGD
  if (length(grep("WGD", unique(tdata$ID))) == 0) {
    including_wgd = FALSE
  }

  tdata$plotpos = -log10(tdata$value)
  tdata$plotpos = (tdata$plotpos - min(tdata$plotpos)) + 0.11

  if (including_wgd) {
    #WGD_linepos <- tdata %>% filter(str_detect(ID,'WGD')) %>% select(ID,value) %>% group_by(ID) %>% summarise(pos=-log10(median(value))) %>% pull(pos)
    WGD_linepos <- tdata %>% filter(str_detect(ID,'WGD')) %>% select(ID,value,plotpos) %>% group_by(ID) %>% summarise(pos=median(plotpos)) %>% pull(pos)
  } else {
    WGD_linepos <- -log10(median(event_medians))
  }
  maxpos <- max(tdata$plotpos)

  midpoint <- median(tdata$plotpos)

  tdata2 <- tdata %>%
    group_by(ID) %>% 
    summarise(mvalue=median(plotpos),max=(max(plotpos)),min=(min(plotpos))) %>% 
    mutate(plotpos=if_else(mvalue>midpoint,min,max),cna='LOH') %>% 
    mutate(ID2=str_remove(ID,' \\[.*')) %>% 
    mutate(ID2=str_remove(ID2,'^.*: '))

  tdata2.1 <- tdata2 %>% filter(mvalue>midpoint)
  tdata2.2 <- tdata2 %>% filter(mvalue<=midpoint)
  tdata2.3 <- tdata %>% group_by(ID,freq) %>% dplyr::slice(1) %>% ungroup()
  tdata2.3$freq <- tdata2.3$freq / 1000

  #event_type_colours = c(`Mut`="#30C358",`Gain`="#DE3E67",`LOH`="#2D4FA5",`HD`="#6FDDEF",`WGD`="#FFD993", `MRCA`="#555555")
  #event_type_colours = c(`Mut`="#30C358",`Gain`="#E60751",`LOH`="#398BE3",`HD`="#664499",`WGD`="#FFA500", `MRCA`="#555555")
  event_type_colours = c(`Mut`="#DD8AB9",`Gain`="#68BFA3",`LOH`="#8E9FCA",`HD`="#F18A62",`WGD`="#664499", `MRCA`="#555555")
  #event_type_colours = c(`Mut`="#FFD993",`Gain`="#DE3E67",`LOH`="#2D4FA5",`HD`="#6FDDEF",`WGD`="#", `MRCA`="#555555")
  event_type_colours = event_type_colours[names(event_type_colours) %in% tdata$cna]
  event_type_labels = c(`Mut`="Mutation",`Gain`="Gain",`LOH`="Loss",`HD`="Homozygous Deletion",`WGD`="WGD",`MRCA`="MRCA")
  event_type_labels = event_type_labels[names(event_type_labels) %in% tdata$cna]

  if(horizontal) {
    gt_angle = 90
    gt_vjust = 0.5
    save_height = 7
    save_width = 10
  } else {
    gt_angle = 0
    gt_vjust=0.5
    save_height = 80/25.3
    save_width = 140/25.3
  }

  text_size <- 1.5

  p2 <- tdata %>% 
    ggplot(aes(x = ID, y = plotpos, fill=cna)) +
    geom_stripes(odd = "#11111106", even = "#00000000")+
    geom_hline(yintercept=WGD_linepos, linetype='dashed', col = ifelse(including_wgd, 'gray70', rgb(1,1,1,0)),linewidth=0.25) +
    geom_violin(color=NA) + #'black',size=0.001
    geom_hline(yintercept=0, linetype='solid', col = 'gray50',linewidth=0.1) +
    geom_hline(yintercept=0.05, linetype='dotted', col = 'gray50',linewidth=0.1) +
    geom_hline(yintercept=0.1, linetype='dotted', col = 'gray25',linewidth=0.1) +
    geom_bar(data=tdata2.3, aes(y=freq), width=2/3, color='black', linewidth = 0.001, stat="identity")+
    geom_text(data=tdata2.3,aes(y=-0.02, label=paste0(format(round(freq*1000,1), nsmall=1),"%")),hjust=1,angle=gt_angle,size=text_size,family = 'Roboto Condensed' )+
    #geom_text(data=tdata2.3,aes(y=0, label=paste0(format(round(freq*1000,1), nsmall=1),"%")),hjust=1,angle=gt_angle,size=text_size,family = 'Roboto Condensed' )+
    stat_summary(fun = "median", geom = "crossbar", width = 0.4,linewidth=0.18, colour = "black", show.legend=FALSE) +
    geom_text(data = tdata2.1,aes(label=ID2),vjust = 0.5,hjust=1,nudge_y = -0.05,size=text_size,angle=gt_angle,family = 'Roboto Condensed' )+
    geom_text(data = tdata2.2,aes(label=ID2),vjust = 0.5,hjust=0,nudge_y = 0.05,size=text_size,angle=gt_angle,family = 'Roboto Condensed' )+
    #labs(x="",y="Ordering scale",fill="Event",size='Frequency (%)') +
    labs(x="",y="",fill="Event",size='Frequency (%)') +
    scale_fill_manual(breaks=names(event_type_colours), values=event_type_colours, labels=event_type_labels) +
    scale_size_binned(n.breaks = 8)+
    guides(fill = guide_legend(byrow = TRUE,override.aes = list(size=2.5,pch=21)),size=guide_legend(override.aes = list(fill='gray50',col='black'))) +
    theme_ipsum_rc(base_family="Roboto Condensed",ticks = FALSE,grid = FALSE,axis_title_just = 'm',axis_title_size = 9)+
    theme(axis.text.y = element_blank()) +
    theme(legend.position="none") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "null"))  +
    scale_y_discrete(expand = c(0.1, 0.1))

  dim1 = 6
  dim2 = 13
  if(horizontal) {
    wid = dim2
    hei = dim1
  } else {
    hei = dim2
    wid = dim1
    p2 <- p2 + coord_flip(clip = 'off')
  }

  #p2

  #ggsave(file=paste0(output_dir, '/', tumour_set, "_PLPlot_test.pdf") ,plot = p2,width = 8,height = 10,device = postscript())
  #ggsave(file=paste0(output_dir, '/', tumour_set, "_PLPlot_test.png") ,plot = p2,width = 8,height = 10,device = png())
  saveRDS(
    tdata %>% dplyr::select(value, ID, plotpos, cna, freq),
    file.path("outputs/02_trajectories", paste0("Fig2", c("b", "c", "d")[trajectory], "_source_data.rds"))
  )
  cairo_pdf(paste0(figure_dir, '/Fig2_', trajectory, "_PLPlotSlim.pdf"), width = wid,height = hei, family="Roboto Condensed")
  print(p2)
  dev.off()
  # ggsave(file=paste0(figure_dir, '/', trajectory, "_PLPlotSlim.pdf") ,plot = p2,width = wid,height = hei,device = cairo_pdf())
  ggsave(file=paste0(figure_dir, '/Fig2_', trajectory, "_PLPlotSlim.png") ,plot = p2,width = wid,height = hei,device = png())

}



