plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
}

bc_rank_plot <- function(stats, raw_cells, filt_cells, save){
  cells <- raw_cells %in% filt_cells
  keep <- !duplicated(stats$total)
  plot_df <- data.frame(Rx=stats$rank, Tx=stats$total, cell=cells)
  plot_df <- plot_df[keep, ]
  png(save, width = 700, height = 500)
  print({
    ggplot(subset(plot_df, plot_df$UMI_Counts >0), aes(x=Barcodes, y=UMI_Counts, col = Cell))+ 
      geom_point(size=3) + 
      geom_hline(yintercept = stats@metadata$knee, lty = 2, col = '#030202', size=1) +
      annotate("text", x=max(plot_df$UMI_Counts), y=stats@metadata$knee+10000, label=NULL, color = "#030202", size=5) +
      geom_hline(yintercept = stats@metadata$inflection, lty = 2, col = '#030202', size = 1) +
      annotate("text", x=max(plot_df$UMI_Counts), y=stats@metadata$inflection+100, label=NULL, color = "#030202", size=5) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
      scale_y_log10(breaks = trans_breaks('log10', function(x) floor(10^x)))+
      labs(x = 'Barcodes', y = 'UMI counts', title ="Kb Count SubLib2 ATCC") +
      guides(colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
      theme_linedraw()+
      theme(legend.title = NULL,
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.text = element_text(size=10),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = c(0.15,0.15))
  })
  dev.off()
}

print_HTML <- function(seq_stats, cell_stats, dir, sample_id){
  system(paste0('base64 ', dir, '/barcode_rank.png > ', dir, '/barcode_rank.txt'))
  b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)
  target <- HTMLInitFile(dir, filename=paste0(sample_id, '_summary'))
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
  HTML("<div class='title'>", file=target)
  HTML.title(' Pre-Processing Summary', HR=1, file = target)
  HTML("</div>", file = target)
  HTML.title(sample_id, HR=2, file = target)
  HTML("<div id='wrapper'>", file=target)
  HTML("<div class='boxed' id='left' align='center'>", file=target)
  HTML.title('Sequencing/Alignment Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', seq_stats$stat, '</td> <td align="right">', seq_stats$value, '</td> </tr>'), file=target)
  HTML('</table> <hr>', file=target)
  HTML.title('Cell Stats', HR=3, file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', cell_stats$stat, '</td> <td align="right">', cell_stats$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML("</div>", file = target)
  HTML("<div class='boxed' id='right' align='center'>", file=target)
  HTML(paste0("<img src='data:', b64_bc, '>'"), file=target)
  HTML("</div>", file = target)
  HTML("</div>", file = target)
  HTML('<style type="text/css">
		.title {
    			background-color: #0972D5;
    			padding: 8px;
    			color: white;
    			position: fixed;
    			top: 0;
    			left: 0;
    			z-index: 999;
    			width: 100%;
		}
		.boxed {
  			border: 1px solid #868D96;
  			padding: 10px;
  			margin: 20px;
		}
		h1 {
			font-family: "Roboto";
			font-size: 33px;
		}
		h2 {
			font-family: "Roboto";
			font-size: 26px;
		}
		h3 {
			font-family: "Roboto";
			font-size: 18px;
		}
		#wrapper {
  			display: flex;
		}
		#left {
  			width: 50%;
		}
		#right {
  			width: 50%;
		}
		table tr:nth-child(even) {
  			background-color: #eee;
		}
		table tr:nth-child(odd) {
  			background-color: #fff;
		}
		table {
  			font-family: "Roboto";
			font-size: 20px;
			border: 1px solid #868D96;
		}
		#mathplayer{
  			height: 80px;
		}
		</style> </head>', file=target)
  HTMLEndFile()
}
