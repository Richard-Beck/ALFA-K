##Some of the data we have looked wrong and was excluded. 
## This function removes those elements from the metadata, so they don't show in the
## already cluttered network plots.  
## The premise is that these files are absent from the processed alfak_inputs folder
prune_children <- function(meta,dir="data/processed/salehi/alfak_inputs"){
  filename_patterns <- paste(meta$datasetname,meta$timepoint,sep="_")
  ff <- list.files(dir)
  c1 <- sapply(filename_patterns,function(fi) sum(grepl(fi,ff))>1) 
  ## the other reason these file patterns can be absent from the data is if they're
  ## the initial "root" passages. All these should be kept.
  ## in general just keeping all passages that have a parent is a sound choice:
  c2 <- sapply(meta$uid,function(idi) idi%in%meta$parent)
  meta[c1|c2,]
}

get_sample_lineage_map <- function(){
  sample_lineage_map <- c(
    SA1035T                   = "SA1035",
    SA1035U                   = "SA1035",
    SA532                     = "SA532",
    SA609                     = "SA609",
    SA609UnBU                 = "SA609",
    SA906a                    = "p53 k.o",
    SA906b                    = "p53 k.o",
    SA039U                    = "p53 w.t",
    SA535_CISPLATIN_CombinedT = "SA535",
    SA535_CISPLATIN_CombinedU = "SA535",
    SA000 = "SA609",
    SA609R2 = "SA609",
    SA535_CISPLATIN_CombinedH = "SA535",
    SA001="SA609",
    SA609H2="SA609",
    SA609H1="SA609",
    SA609R2b="SA609",
    SA609TIIDT="SA609"
  )
}

## maps a lineage name (fi) to a cell line.
## lineage assumed to be in format sampleName_timepoint_l_other_details
assign_labels <- function(fi,meta){
  sample_lineage_map <- get_sample_lineage_map()
  lut <- meta$datasetname
  names(lut) <-  paste(meta$datasetname,meta$timepoint,sep="_")
  id <- strsplit(fi,split="_l_") |> unlist() |> head(1)
  id_intermediate <- lut[id]
  sample_lineage_map[id_intermediate]
}

## df <- read.csv("data/raw/salehi/metadata.csv") 
basic_lineage_plot <- function(df,  plot_pars = list(base_node_color="black",
                                                      base_edge_color="black")){
  df <- prune_children(df)
  
  # Create dummy root to fix missing parent values
  dummy_uid <- "ROOT"
  dummy_row <- df[1, ]
  dummy_row$uid <- dummy_uid
  dummy_row$datasetname <- "ROOT"
  df_fixed <- rbind(df, dummy_row)
  df_fixed$parent[is.na(df_fixed$parent) | df_fixed$parent == ""] <- dummy_uid
  
  df_fixed <- cbind(df_fixed[, c("uid", "parent")],
                    df_fixed[, !colnames(df_fixed) %in% c("uid", "parent")])
  
  # Build edge list and graph ----------------------------------------------
  edges_df <- df_fixed %>% 
    filter(uid != dummy_uid) %>% 
    mutate(from_dummy = (parent == dummy_uid)) %>% 
    select(from = parent, to = uid,  from_dummy)
  
  g <- graph_from_data_frame(d = edges_df, vertices = df_fixed, directed = TRUE)
  
  tg <- as_tbl_graph(g) %>% 
    mutate(dummy = (name == dummy_uid),
           node_type = ifelse(dummy, "dummy", ifelse(parent == dummy_uid, "root", "non-root")))
  
  # Define layout and compute label positions -------------------------------
  set.seed(42)
  layout <- create_layout(tg, layout = "dendrogram", circular = TRUE,
                          height = -node_distance_to(1, mode = "all"))
  
  p_network <- ggraph(layout) + 
    geom_edge_link(aes(linetype=ifelse(from_dummy, "blank", "solid")),
                   edge_width = 0.5,show.legend = F,color=plot_pars$base_edge_color) +
    # Plot dummy nodes as white points
    geom_node_point(data = filter(layout, dummy),
                    aes(x = x, y = y), color = "white") +
    geom_node_point(data = filter(layout, !dummy),
                    aes(x = x, y = y), color=plot_pars$base_node_color) +
    scale_edge_linetype_manual(values = c(solid = "solid", blank = "blank"),
                               guide = "none")+
    theme_classic()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA)
    )
  
  
  list(plt=p_network,layout=layout)
}


generate_segment_highlights <- function(plt,lins){
  
  
  highlights <- do.call(rbind,lapply(1:length(lins),function(i){
    uid <- sapply(strsplit(lins[[i]],"-"),tail,1)
    parent <- sapply(strsplit(lins[[i]],"-"),head,1)
    pdx <- assign_labels(names(lins)[i],meta = m)
    data.frame(name=uid,parent,pdx,traj=i,row.names = NULL)
  }))
  highlights <- split(highlights,f=highlights$pdx)
  highlights <- do.call(rbind,lapply(highlights,function(hi){
    hi$categorical_variable <- LETTERS[c(1,1+cumsum(diff(hi$traj)>0))]
    hi
  }))
  # First: Join highlights to layout to get edge coordinates
  # layout$name gives node id; layout$parent is the parent id
  
  edge_highlights <- highlights %>%
    rename(to = name) %>%
    left_join(plt$layout %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x_to = x, y_to = y) %>%
    left_join(plt$layout %>% select(name, x, y), by = c("parent" = "name")) %>%
    rename(x_from = x, y_from = y)
}