###############################################################################
##########################################
####Project resistance in P.Aeruginosa####
 ##########################################
###############################################################################

### Carga de librerias ###
devtools::install_github("rstudio/addinexamples", type = "source")
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(rvest)
library(stringr)
library(imputeTS)
library(factoextra)
library(FactoMineR)
library(GGally)
library(ggpubr)
library(ggplot2)
library(ggfortify)
library(ggcorrplot)
library(dslabs)
library(tseries)
library(ggpubr)


### Loading and bulding datasets for first part (R_index)
Pseudomona_db <- read.csv("data_complete_int.csv", header=TRUE) %>% 
    rename('R_multi' = 'multi_R_perc', 
           'R_index' = 'r_index')

### creating a function that become values above 0 to 0 (because this are frequencies)
min_0 <- function(x){ifelse(x < 0, 0, x)}

### applying function min_0 to every importan resistance variable ###
Pseudomona_db <- Pseudomona_db %>% mutate(across(c("Aminoglycosides", "Carbapenems", "Ceftazidime",
                                    "Fluoroquinolones", "Piperacilina_taz", "R_multi"),
                                  min_0))

### creating chunks of data to further analysis & generating data variables ###
Pseudomona <- Pseudomona_db[, -c(3, 11:12, 14:26)]

p.aeuro_r <- Pseudomona

p.aeuro_r_2 <- cbind(p.aeuro_r, Pseudomona_db$R_index/100) %>% rename("R_index" ="Pseudomona_db$R_index/100")

Data <- as_tibble(Pseudomona[, -c(1:4, 11:12, 14:26)])

n <- nrow(Data) 

### Graphic analysis of resistance variables
gath_p.aeuro <- Pseudomona %>% 
    gather(key = "antibiotic", value = "resistance", 5:10)

### Density plot all resistance variables 
gath_p.aeuro %>% ggplot(aes(resistance)) + 
    scale_x_continuous(trans = "log2",
                       labels = scales::number_format(accuracy = 0.01)) +
    geom_density(aes(fill=factor(antibiotic)), 
                 alpha=0.5, 
                 position = "stack") + 
    labs(title="Density plot", 
         subtitle="Antibiotic resistance total",
         #caption="Source: mpg",
         x="Resistance Log2(%)",
         fill="Antibiotics") + 
    ds_theme_set(legend.title = element_blank(), legend.position="bottom")

### Box plot 
### BP by Region 
gath_p.aeuro %>%
    mutate(Region = reorder(Region, resistance, FUN = mean)) %>%      # reorder
    ggplot(aes(Region, resistance, fill = antibiotic)) +    # color by continent
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Resistance Log2(%)") +
    xlab("") + 
    geom_jitter(width = 0.1, alpha = 0.1)+
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01)) + 
    ggtitle("Boxplot by Region")+ 
    ds_theme_set(legend.title = element_blank(), legend.position="bottom")

### BP by Country
gath_p.aeuro %>%
    mutate(Country = reorder(Country, resistance, FUN = mean)) %>%      # reorder
    ggplot(aes(Country, resistance, fill = antibiotic)) +    # color by antibiotic
    geom_boxplot() +
    ylab("Resistance Log2(%)")+
    xlab("") + #geom_point(show.legend = FALSE) + 
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01)) + 
    ggtitle("Boxplot by country")+ 
    ds_theme_set(legend.title = element_blank(), legend.position="bottom", 
                 axis.text.x = element_text(angle = 90, hjust = 1))

### Timeserie all resistance variables
Serie <- p.aeuro_r %>% group_by(Year) %>% summarise_if(is.numeric, mean)
Serie_index <- p.aeuro_r_2 %>% 
    group_by(Year) %>% 
    select(R_multi, R_index) %>% 
    summarise_if(is.numeric, mean)

tsbact <- ts(Serie[-1],
             frequency = 1, 
             start = 2005)

tsbact_2 <- ts(Serie_index[-1],
               frequency = 1, 
               start = 2005)

### Graph time serie
autoplot(tsbact, facets=FALSE, label = TRUE) +
    geom_line(size= 1)+
    labs(title="Time Series plot ", 
         subtitle="Antibiotics & multiresistance", 
         caption="Source: Prepared by the authors", 
         y="Resistance (%)", 
         color=NULL) +
    #   geom_text(data = labels, aes(x, y, label = names(serie[-1])), size = 3) +
    #   theme(legend.position = "none")
    #aes(linetype = plot_group,
    #   size = plot_group) +  
    scale_size_manual(values = c(1, 1, 1, 1.5, 1, 1)) + 
    scale_x_continuous(labels = p.aeuro_r$Year, breaks = p.aeuro_r$Year)+# title and caption
    ds_theme_set(plot.caption = element_text(),
                 legend.title = element_blank(), legend.position="bottom") #este theme solicita dslabs


### Graph time serie R_index Vs R_multi ###
autoplot(tsbact_2, facets=FALSE, label = TRUE) +
    geom_line(size= 1)+
    labs(title="Time Series plot ", 
         subtitle="Resistance index & multiresistance", 
         caption="Source: Prepared by the authors", 
         y="Resistance (%)", 
         color=NULL) +
    #   geom_text(data = labels, aes(x, y, label = names(serie[-1])), size = 3) +
    #   theme(legend.position = "none")
    #aes(linetype = plot_group,
    #   size = plot_group) +  
    scale_size_manual(values = c(1, 1, 1, 1.5, 1, 1)) + 
    scale_x_continuous(labels = p.aeuro_r$Year, breaks = p.aeuro_r$Year)+# title and caption
    ds_theme_set(plot.caption = element_text(),
                 legend.title = element_blank(), legend.position="bottom") #este theme solicita dslabs

Corr_resist <- round(cor(Data), 2)

### Correlation graph resistance variables
ggcorrplot(Corr_resist, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"),
           ggtheme=theme_bw, title = "Correlogran of antibiotics resistances")

#####################################################################
### comparison data between R_index and R_multi by visualization ####
#####################################################################

graph_multi <- p.aeuro_r_2 %>% select(Bacteria, Year, Country, Region, R_multi, R_index) %>% 
    gather(key = "antibiotic", value = "resistance", 5:ncol(.))

### BP by Region R_index and R_multi
graph_multi %>%
    mutate(Region = reorder(Region, resistance, FUN = mean)) %>%      # reorder
    ggplot(aes(Region, resistance, fill = antibiotic)) +    # color by continent
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Resistance Log2(%)") +
    xlab("") + 
    geom_jitter(width = 0.1, alpha = 0.1)+
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01)) + 
    ggtitle("Boxplot by Region")+ 
    ds_theme_set(legend.title = element_blank(), legend.position="bottom")

### BP by Country R_index and R_multi
graph_multi %>%
    mutate(Country = reorder(Country, resistance, FUN = mean)) %>%      # reorder
    ggplot(aes(Country, resistance, fill = antibiotic)) +    # color by antibiotic
    geom_boxplot() +
    ylab("Resistance Log2(%)")+
    xlab("") + #geom_point(show.legend = FALSE) + 
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01)) + 
    ggtitle("Boxplot by country")+ 
    ds_theme_set(legend.title = element_blank(), legend.position="bottom", 
                 axis.text.x = element_text(angle = 90, hjust = 1))

###########################################
### K-means method for every antibiotic ### 
###########################################

mojena = function(hc){
    n_hd = length(hc$height)
    alp_g = 0 ; alpha = hc$height[n_hd:1]
    for(i in 1:(n_hd-1)){
        alp_g[i] = mean(alpha[(n_hd-i+1):1])+1.25*sd(alpha[(n_hd-i+1):1])
        # alp_g[i] = mean(alpha[(n_hd-i+1):1])+3.5*sd(alpha[(n_hd-i+1):1])
    }
    nog = sum(alp_g<= alpha[-n_hd]) + 1
    plot(alpha[-n_hd], pch=20, col=(alp_g>alpha[-n_hd])+1, main = paste("Número óptimo de grupos =",nog),
         ylab = expression(alpha[g]), xlab="Nodos")}

attach(Pseudomona_db)
mean_pseudo<- Pseudomona_db %>% group_by(Country) %>%
    summarise(Amino = mean(Aminoglycosides), 
              Carba = mean(Carbapenems), 
              Fluo = mean(Fluoroquinolones), 
              Cefta = mean(Ceftazidime), 
              Piper = mean(Piperacilina_taz))

multi_pseudo <- as_tibble(mean_pseudo[, -1])
rownames(multi_pseudo) <- mean_pseudo$Country

distance<- get_dist(multi_pseudo, stand = TRUE, method = "euclidean") 
fviz_dist(distance, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#Dendograma método de average
average_pseudo= hclust(distance, method = "average")
fviz_dend(average_pseudo, k = 3, # Cortar en 3 grupos
          cex = 0.5, # tama?o de las etiquetas
          k_colors = "npg", # color de los grupos
          color_labels_by_k = TRUE, # color de las etiquetas de los grupos
          rect = TRUE # Agregar un rect?ngulo a los grupos
)

#hacemos los cluster
cluster1 <- scale(multi_pseudo)
set.seed(123) ## Semilla
km.res_p.aeruginosa <- kmeans(cluster1, 3,nstart = 25)
cluster2 <- cbind(multi_pseudo, 
                  cluster1=km.res_p.aeruginosa$cluster)

fviz_nbclust(cluster1,kmeans,method = "wss")+
    geom_vline(xintercept = 3, linetype=2) ## Advertencia

## Visualizaci?n conglomerados K-medias
fviz_cluster(km.res_p.aeruginosa, 
             data = cluster1, 
             palette="jco",
             labelsize = 8,
             star.plot= TRUE,
             pointsize = 2,
             repel = TRUE,
             ggtheme=theme_minimal()
)


################################################################################
### loading data set for PA resistance, resistance index, and other variables### 
################################################################################

data_PA <- read.csv("data_complete_intGOV.csv")

### setting and recalling a new dataset including only interes variables ###
vars_corr <- data_PA %>% select(r_index, multi_R_perc, DDD_sys_commun, 
                                per_cap_US., GDP_total, GDP_._health, 
                                Out_pocket_exp, HDR, ctrl_corrup, GOV_effect, 
                                rule_law, rural_pop) %>% 
    rename('R_index' = 'r_index',
           'R_multi' = 'multi_R_perc',
           'GDP_Health' = 'GDP_._health', 
           'HDI' = 'HDR')

### function of correlation
corr_vars <- round(cor(vars_corr), 2)

### Correlation graph all variables
ggcorrplot(corr_vars, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 2, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlogram of variables", 
           ggtheme=theme_bw)



Data_index <- read.csv("data_complete_intGOV.csv")
Data_index<- Data_index %>% mutate(r_index = r_index/100) 

Data_index %>% write.csv(file = 'data_index_GOV.csv', row.names = F) ## creación de archivo .csv con la BD interpolada

#### Graph effects country & year, R_multi ###

eff_country <- read.csv2("country_effect_Rmulti.csv") %>% set_colnames(c('Country', 'Effect'))


eff_country_Lv <- eff_country %>% mutate(Effect = as.numeric(as.character(Effect))) %>% 
    mutate(Effect = round(Effect, 2)) %>% 
    mutate(., Effects_level = ifelse(Effect <=0, "Low", "High")) %>% arrange(Effect)

eff_country_Lv$`Country` <- factor(eff_country_Lv$`Country`, levels = eff_country_Lv$`Country`) 
eff_country_Lv$`Effects_level` <- factor(eff_country_Lv$`Effects_level`) 

### DOT PLOT using GGPLOT
ggplot(eff_country_Lv, aes(x=`Country`, y=Effect, label=Effect)) + 
    geom_point(stat='identity', aes(col=Effects_level), size=8)  +
    scale_color_manual(name="Effects_level", 
                       labels = c("High", "Low"), 
                       values = c("High"="Blue", "Low"="Red")) + 
    geom_text(color="white", size=3, fontface = "bold") +
    labs(title="Diverging Dot Plot", 
         subtitle="Normalized mileage from 'mtcars': Dotplot") + 
    ylim(-0.25, 0.25) +
    coord_flip() + 
    theme_minimal()

### DOT PLOT using ggdotchart
require(ggpubr)
ggdot_rmulti <- ggdotchart(eff_country_Lv, x = "Country", y = "Effect",
           color = "Effects_level",                                # Color by groups
           palette = c("#00AFBB","#FC4E07"), # Custom color palett                      # Sort value in descending order
           add = "segments",                           # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 1.5), # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "Effects_level",                                # Order by groups
           dot.size = 8,                                 # Large dot size
           label = (eff_country_Lv$Effect),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 7,
                             vjust = 0.5),             # Adjust label parameters
           title = "Dot Plot R_multi by Country",
           subtitle = "Effects per Country in data panel model R_multi",
           caption = "Source: Prepared by the authors", 
           ggtheme = theme_pubr()                        # ggplot2 theme
)

ggdot_rmulti <- ggdot_rmulti +
    font("xlab", size = 12)+
    font("ylab", size = 12)+
    font("xy.text", size = 10)

#### Graph effects country, R_index ###

eff_index <- read.csv2("country_effect_Rindex.csv") %>% set_colnames(c('Country', 'Effect'))


eff_index_Lv <- eff_index %>% mutate(Effect = as.numeric(as.character(Effect))) %>% 
    mutate(Effect = round(Effect, 2)) %>% 
    mutate(., Effects_level = ifelse(Effect <0, "Low", "High")) %>% arrange(Effect)

eff_index_Lv$`Country` <- factor(eff_country_Lv$`Country`, levels = eff_country_Lv$`Country`) 
eff_index_Lv$`Effects_level` <- factor(eff_country_Lv$`Effects_level`) 

### DOT PLOT using GGPLOT
ggplot(eff_index_Lv, aes(x=`Country`, y=Effect, label=Effect)) + 
    geom_point(stat='identity', aes(col=Effects_level), size=8)  +
    scale_color_manual(name="Effects_level", 
                       labels = c("High", "Low"), 
                       values = c("High"="Blue", "Low"="Red")) + 
    geom_text(color="white", size=3, fontface = "bold") +
    labs(title="Diverging Dot Plot", 
         subtitle="Normalized mileage from 'mtcars': Dotplot") + 
    ylim(-25, 45) +
    coord_flip() + 
    theme_minimal()

### DOT PLOT using ggdotchart
require(ggpubr)
ggdot_rindex <- ggdotchart(eff_index_Lv, x = "Country", y = "Effect",
           color = "Effects_level",                                # Color by groups
           palette = c("#00AFBB","#FC4E07"), # Custom color palett
           #sorting = "ascending",                        # Sort value in descending order
           add = "segments",                           # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 1.5), # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "Effects_level",                                # Order by groups
           dot.size = 8,                                 # Large dot size
           label = (eff_index_Lv$Effect),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 7,
                             vjust = 0.5),# Adjust label parameters
           title = "Dot Plot R_index by Country",
           subtitle = "Effects per Country in data panel model R_index",
           caption = "Source: Prepared by the authors",
           ggtheme = theme_pubr()# ggplot2 theme
)

ggdot_rindex <- ggdot_rindex +
    font("xlab", size = 12)+
    font("ylab", size = 12)+
    font("xy.text", size = 10)


