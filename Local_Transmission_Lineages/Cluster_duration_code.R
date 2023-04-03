## sliding mean cluster length plot by region
require(cowplot)
require(ggplot2)
require(tidyr)
require(cowplot)
require(dplyr)
require(tibble)
require(stringr)
require(lubridate)
require(janitor)

# takes downstreamtree file from the LTL script

alldays <- seq(ymd(rootdate),ymd(mrsd), by = "days")

# start from first sample
earliest_date <- as.Date(min(date_decimal(sapply( strsplit( tre$tip.label, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )}))))
alldays <- seq(earliest_date,ymd(mrsd), by = "days")


# Time to sampling ---------------------------------------------------------

detection <- do.call(rbind.data.frame, lapply(downstreamtree$V8, function(x) min(as.numeric(gsub("_", "",str_extract_all(x, "\\_\\d+\\.?\\d*", simplify = T)))))) ## more complex as some dates don't have decimals
names(detection) <- "detection"
downstreamtree$detection <- as.Date(date_decimal(detection$detection))
downstreamtree$ttd <- downstreamtree$detection - downstreamtree$V3

## Mean time to detection ----
rolling_mean_ttd <- function(date_end, downstreamtree, roll_period){
  order_regi <- levels(as.factor(downstreamtree$V2))
  day_include <- seq(date_end-roll_period, date_end, by = "days")
  suppressMessages({
    matchdata <-  downstreamtree  %>%
      drop_na() %>%
      rowwise() %>%
      mutate(match = ifelse(any(between(day_include, V3, V4)), 1, 0)) %>% # based on when extant rather than time to first detection
      filter(match == 1)
    
    matchdata <- matchdata %>%
      group_by(V2)%>%
      summarise(mean = mean(ttd, rm.na=T),
                n = n(),
                s = sd(ttd),
                error = 1.96*(s/sqrt(n)),
                lower = ifelse(mean-error>=0,mean-error,0),
                upper = mean+error) %>%      
      right_join(expand.grid(V2 = order_regi), by = "V2") %>%
        mutate(date = date_end)
  })
  
  return(matchdata[order(match(names(matchdata), order_regi))])
}
rolling_mean_ttd(as.Date("2022-07-01"), downstreamtree, roll_period = 730)

mean_ttd <- lapply(alldays,rolling_mean_ttd, downstreamtree=downstreamtree, roll_period = 1)

ttd_plot <- do.call(rbind.data.frame, mean_ttd)

ttd_mean<- ggplot() +
  geom_line(data = ttd_plot, aes(x = date, y = mean, color = V2), show.legend = F)+
  geom_ribbon(data = ttd_plot, aes(ymax = upper, ymin = lower, x = date, fill = V2), alpha = 0.5, show.legend = F)+
  scale_color_manual(values = getPalette, drop=FALSE, na.value="white")+
  scale_fill_manual(values = getPalette, drop=FALSE, na.value="white")+
  geom_hline(yintercept = 547, color = "black")+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(0,700))+
  theme_bw()+
  scale_x_date(limits = c(as.Date("2014-01-01"), ymd(mrsd)), expand = c(0,0))+
  labs(x = "Date",
       y = "Rolling mean time to detection (days)") +
  facet_wrap(~V2, ncol = 2)

ttd_mean



# Time to sampling --------------------------------------------------------


downstreamtree$onwards <- do.call(rbind.data.frame, lapply(downstreamtree$V8, function(x) str_detect(x, "collapse")))

downstreamtree$onwards <- downstreamtree$onwards[,1]

downstreamtree <- downstreamtree %>%
  mutate(Category = case_when(diff<=730 & V6<5 & onwards == F ~ "dead end",
                              diff>=365 & V6>10 ~ "persistent",
                              onwards == T ~ "export",
                              T ~ "other")) 

downstreamtree %>% group_by(V2) %>% tabyl(Category)


table(downstreamtree$V2, downstreamtree$Category)


rolling_type_clusters <- function(date_end, downstreamtree, roll_period){
  order_regi <- levels(as.factor(downstreamtree$V2))
  day_include <- seq(date_end-roll_period, date_end, by = "days")
  suppressMessages({
    matchdata <-  downstreamtree  %>%
      drop_na() %>%
      rowwise() %>%
      mutate(match = ifelse(any(between(day_include, V3, V4)), 1, 0)) %>%
      filter(match == 1)
    
    matchdata <- matchdata %>%
      group_by(V2, Category) %>%
      tally() %>% 
      mutate(freq = prop.table(n)) %>%
      right_join(expand.grid(V2 = order_regi), by = "V2") %>%
      mutate(date = date_end)
  })
  
  return(matchdata[order(match(names(matchdata), order_regi))])
}

rolling_type_clusters(as.Date("2022-03-14"), downstreamtree, roll_period = 1)


cluster_types <- lapply(alldays,rolling_type_clusters, downstreamtree=downstreamtree, roll_period = 1)

count_cluster_types <- do.call(rbind, cluster_types)
                
count_cluster_types$Category <- relevel(as.factor(count_cluster_types$Category), 'other')  



ggplot() +
  geom_bar(data = count_cluster_types, aes(x = date, y = n, fill = Category, color = Category), stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("gray80", "orange", "deeppink", "purple4"))+
  scale_fill_manual(values = c("gray80", "orange", "deeppink", "purple4"))+
  theme_bw()+
  scale_x_date(limits = c(earliest_date, ymd(mrsd)), expand = c(0,0))+
  labs(x = "Date",
       y = "Number of chains") +
  facet_wrap(~V2, ncol = 2) +
  theme(legend.direction = "horizontal", legend.position = "bottom") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5), color = guide_legend(title.position = "top", title.hjust = 0.5)) 


ggsave("chain_types_noother.png", width = 8, height = 8, dpi = 700)



