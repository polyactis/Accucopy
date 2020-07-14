#!/usr/bin/env Rscript
args=commandArgs();
cat("Arguments are", args, "\n")

infer_out_df=readLines(args[6])
infer_details_df=readLines(args[7])
autocor_df=read.delim(args[8],head=TRUE)
output_filepath=args[9]

library(ggplot2)
library(grid)
library(scales)

new_theme <- theme_set(theme_bw())
new_theme <- theme_update(
    plot.title = element_text(colour="black",size=28,angle=0,hjust=0.5,vjust=1,face="plain"),
    plot.margin = unit(c(.25, .25, .25, .25),'in'), # top, right, bottom, left
    axis.title.x = element_text(colour="black",size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
    axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="plain")
)

theme_set(theme_bw())
theme_set(new_theme)

y_axis_value_for_text=seq(1,0,-0.06)
result_summary_text=c(infer_out_df,"","details",infer_details_df)
result_summary_text = gsub("\t", "  ", result_summary_text)
t_len=length(result_summary_text)


#pdf(file=output_filepath, onefile=TRUE, paper="special",width=24,height=20)
jpeg(file=output_filepath, quality=100, width=1024, height=1024)

vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
grid.newpage()
pushViewport(viewport(layout=grid.layout(13, 10)))

cat("Plotting summary result ...\n")
tx=data.frame(x=rep(0,t_len),y=y_axis_value_for_text[1:t_len],lb=result_summary_text)
p_summary_text=ggplot()+
    geom_text(aes(x, y, label=lb),data=tx,hjust=0,size=5)+
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(limits=c(0,1))+
    element_blank()

cat("Plotting rc_ratio autocor ...\n")
p_rc_ratio_autocorr=ggplot(autocor_df, aes(read_count_ratio, correlation)) + geom_point() +
    xlab("TRE")+
    ylab("Auto-correlation")

cat("Plotting rc_ratio autocor ...\n")
autocor_pos_df = autocor_df[autocor_df['correlation']>0,]
autocor_pos_df['log_cor'] = log10(autocor_pos_df['correlation'])
p_log10_autocorr=ggplot(autocor_pos_df, aes(read_count_ratio, log_cor)) + geom_point() +
    xlab("TRE")+
    ylab("log10(Auto-correlation)")

print(p_summary_text, vp=vplayout(1:3,1:10))
print(p_rc_ratio_autocorr, vp=vplayout(4:8,1:10))
print(p_log10_autocorr, vp=vplayout(9:13,1:10))

dev.off()
