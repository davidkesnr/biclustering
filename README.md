Biclustering multiple variables with a semi-automated cluster
optimization protocol
================
David Kesner

required packages:

``` r
library(reshape2)
library(parallel)
library(foreach) 
library(doParallel)
library(biclustermd)
library(grDevices)
library(plot.matrix)
library(RColorBrewer)
library(tidyr)
library(rnaturalearth)
```

## Background

In this analysis, I attempt to define a European spatial stratification
that separates regions of coherent environmental change with specific
relevance for wildfire activity over multiple millennia. To do this, I
use a biclustering algorithm to group spatiotemporal environmental
datasets spanning the Holocene epoch. The variables chosen to represent
fire-relevant environmental change are summer temperature and summer
precipitation anomalies, representing fire season dryness, from the
Mauri et al. (2015) dataset, as well as forest cover percentage data
from the Zanon et al. (2018) dataset, representing vegetation available
to burn as fuel. The three datasets have been processed from their
original netcdf format into .csv files at the same spatial and temporal
resolution, and the forest cover values have been anomalized to match
the Mauri et al. (2015) climate data.

Let’s visualize the forest cover data:

``` r
#read in forest data
forest <- read.csv("forest.csv")
head(forest)
```

    ##   ID_ENTITY year       vals  lon  lat
    ## 1         1  100  0.0000000 -0.5 31.5
    ## 2         1 1000  5.5406794 -0.5 31.5
    ## 3         1 2000 -4.9102641 -0.5 31.5
    ## 4         1 3000  0.8845382 -0.5 31.5
    ## 5         1 4000  6.1101104 -0.5 31.5
    ## 6         1 5000  1.9467614 -0.5 31.5

``` r
unique(forest$year)
```

    ##  [1]   100  1000  2000  3000  4000  5000  6000  7000  8000  9000 10000 11000
    ## [13] 12000

The data is structured with multiple time points per grid cell, and each
grid cell has an associated ID in the `ID_ENTITY` column. We can see
that the data is at a 1000-year time step spanning the 0-12000BP
interval and is expressed as anomalies relative to a 100
years-before-present (BP) value. Lets plot one time slice:

``` r
#get the world outline shapefile:
world <- ne_countries(scale = "medium", returnclass = "sf")

#get one time slice
for.slice <- forest[forest$year %in% "3000",]

#plot
ggplot() +
  geom_point(data = for.slice, aes(x = lon, y = lat, fill = vals),shape = 21, size = 6)+
  scale_fill_gradient2(low="red", high="forestgreen", mid = 'white', midpoint = 0)+
  scale_x_continuous(limits = c(-11,44), breaks = seq(-10, 40, by=10)) +
  theme(legend.text=element_text(size=22),legend.position = 'bottom',legend.spacing.x = unit(0.5,'cm'),legend.spacing.y = unit(0.5,'cm'),  legend.key.height = unit(2, 'cm'), legend.key.width = unit(6, 'cm'), legend.title=element_blank() , panel.background = element_rect(fill = "aliceblue"), axis.text = element_text(size=18, face='bold', colour = 'black'), axis.title = element_text(size=22, face='bold', colour = 'black'), title = element_text(size = 10, colour = 'black'))+
  scale_y_continuous(limits = c(29, 73), breaks = seq(30, 90, by=10)) +
  labs(x = "", y =  "")+
  geom_sf(data = world, fill="transparent", color="grey80")
```

![](biclustering_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

We can see that forest cover was generally higher across Europe at
3000BP relative to 100BP. We also learn from this plot that the data
covers a range of roughly 30-70 degrees north in latitude, and has a
longitudinal range from 10 degrees west to about 43 degrees east. The
climate data has the same data coverage.

## Defining the stratification

The goal is to group the grid cells of the three datasets based on the
similarity in their environmental variation over time. To do this, I use
a biclustering algorithm, namely the algorithm developed by Reisner et
al. (2020). This algorithm is convenient for my purposes due to its
ability to handle missing data (unlike other clustering algorithms
e.g. k-means). It also provides the ability to define temporal clusters,
which is useful as it allows for higher precision in the assignment of
grid cells to spatial clusters.

One problem is that the clustering algorithm only allows clustering of a
single variable at a time, but I am attempting to use three variables at
once. To solve this, I have concatenated the three variables into a
single dataframe, aligning them by their grid-cell ID’s. This creates a
dataframe with 13x3 time points.

There are a few problems that this introduces. One is that the variables
are in different units: summer temperature and precipitation are
provided as anomalies relative to the 100BP time bin, expressed in
degrees Celsius and mm/month, and forest cover is provided as
percentages. This makes the three variables span orders of magnitude.
Clustering them in their native units will result in the variables on
higher orders of magnitude to reduce the visibility of the variation in
variables on lower orders of magnitude to the algorithm. To circumvent
this, I have scaled the three variables to have identical interquartile
ranges.

I have also removed the 100BP time points because they contain no
variation and will not influence the clustering results.

lets look at the data:

``` r
#read in concatenated dataset
bcdat <- read.csv("bcdat.csv")
head(bcdat)[1:5]
```

    ##    tsumm_1000 tsumm_2000  tsumm_3000  tsumm_4000 tsumm_5000
    ## 1  0.42269474 -1.6206330  0.06331801  0.41399574 -1.1154922
    ## 2 -0.04767931 -1.2627360 -0.43620682 -0.01457942 -1.0451154
    ## 3 -0.27145228 -1.0315173 -0.56234670 -0.18184218 -1.0622907
    ## 4 -0.47570041 -0.8854572 -0.64877009 -0.30691746 -1.1581845
    ## 5  0.69507378 -1.3439040  0.79554415  0.50073707 -1.0228367
    ## 6  0.68007815 -0.9456235  0.47323859  0.29046571 -0.8679315

``` r
names(bcdat)
```

    ##  [1] "tsumm_1000"   "tsumm_2000"   "tsumm_3000"   "tsumm_4000"   "tsumm_5000"  
    ##  [6] "tsumm_6000"   "tsumm_7000"   "tsumm_8000"   "tsumm_9000"   "tsumm_10000" 
    ## [11] "tsumm_11000"  "tsumm_12000"  "psumm_1000"   "psumm_2000"   "psumm_3000"  
    ## [16] "psumm_4000"   "psumm_5000"   "psumm_6000"   "psumm_7000"   "psumm_8000"  
    ## [21] "psumm_9000"   "psumm_10000"  "psumm_11000"  "psumm_12000"  "forest_1000" 
    ## [26] "forest_2000"  "forest_3000"  "forest_4000"  "forest_5000"  "forest_6000" 
    ## [31] "forest_7000"  "forest_8000"  "forest_9000"  "forest_10000" "forest_11000"
    ## [36] "forest_12000"

We can see that the columns are labelled according to the time point as
well as the variable from which they were taken. The row names match the
grid cell ID’s, which is important for later plotting.

## Biclustering

The biclustering algorithm requires the user to arbitrarily select the
appropriate number of spatial and temporal clusters. However, I would
like to find an optimal combination of clusters based on the variation
inherent to the dataset. To do this, I have created a biclustering
optimization algorithm, which automatically selects an optimal site/time
cluster combination from a range of pre-specified combinations.

This algorithm makes use of an output of the biclustering algorithm
which measures the reduction in the sum of squared errors across the
dataset that a given clustering implementation achieves. This measure
naturally tends towards 0 as the number of parameters estimated
increases, so it does not assist in selecting the most parsimonious
model for the data. This is why I combine this metric with a modified
Aikaike Information Criterion (AIC) formula to create the final
optimization algorithm. The optimization algorithm therefore rewards
parsimony by penalising biclustering implementations in proportion to
the number of clusters used, while rewarding implementations for the
reduction in dataset noise that they achieve. This is the formulation:

![
\\begin{equation}
k = n.ln⁡\\Biggl(\\bigl(\\frac{sse}{n} \\bigr)+2(r.c)\\Biggr)
\\end{equation}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Bequation%7D%0Ak%20%3D%20n.ln%E2%81%A1%5CBiggl%28%5Cbigl%28%5Cfrac%7Bsse%7D%7Bn%7D%20%5Cbigr%29%2B2%28r.c%29%5CBiggr%29%0A%5Cend%7Bequation%7D%0A "
\begin{equation}
k = n.ln⁡\Biggl(\bigl(\frac{sse}{n} \bigr)+2(r.c)\Biggr)
\end{equation}
")

![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is the number of grid cells,
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r "r")
and
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c")
are the number of spatial and temporal clusters respectively, and
![sse](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;sse "sse")
is the final sum of squared errors across all biclusters within a given
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination. The
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination that results in the lowest value of
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
is selected as the optimal combination. We can now use this algorithm to
select from multiple
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combinations, and see which minimizes this value
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k").

Since the biclustering algorithm can be computationally intensive, and I
run it tens or possibly hundreds of times to find the best
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination, this section of the analysis is best done by using a loop,
which runs multiple loop iterations simultaneously by distributing them
across multiple CPU’s. To do this, I use the `foreach` and `parallel`
packages. First, I register the number of cores to recruit for this
purpose. I leave one core for available for general computation outside
of R:

``` r
numCores <- detectCores()-1 #get number of cores
numCores 
```

    ## [1] 15

I am recruiting 15 CPU’s for the biclustering.

The `foreach` loop uses a single index vector as a loop counter, but I
need to be able to loop through two indices (one for
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r "r")
and one for
![c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c "c"))
simultaneously. To achieve this, I make a dataframe with all
combinations of individual
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
pairs, and will use the row index of the dataframe as the `foreach` loop
counter to sequentially extract and run the
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination contained in each row. I set the range of site clusters to
attempt from 2 to 18, and the time clusters from 2 to 13:

``` r
#get an index dataframe
idx <- data.frame(col_t = rep(2:18,each = 12), #col_t denotes time points/columns
                  row_s = rep(seq(2,13,1),17)) #row_s denotes sites/rows
```

The dataframe should have all combinations of 2-18 site and 2-13 time
clusters. Lets check this:

``` r
idx[10:15,]
```

    ##    col_t row_s
    ## 10     2    11
    ## 11     2    12
    ## 12     2    13
    ## 13     3     2
    ## 14     3     3
    ## 15     3     4

``` r
nrow(idx)
```

    ## [1] 204

There are 204 rows, which agrees with all possible combinations (12x17),
and the alignment of the site and time clusters appears correct, with
each site cluster value having 2-13 time clusters.

I am ready to run the `foreach` loop. I’ll calculate the
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
metric for each loop iteration, and for this I need to have the total
number of sites
(![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n"))
on hand:

``` r
nsites <- length(rownames(bcdat))
```

Now, I specify the way that `foreach` must combine the output
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
values of each iteration running in parallel. I choose `cbind` because I
want outputs to be concatenated together into one vector that matches
the order of rows in the `idx` dataframe. This makes it simple to add
the loop output as a row to the `idx` dataframe.

I run the biclustering algorithm with 500 repetitions for each
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination to avoid convergence to local optima (see Reisner et
al. 2020 for more information). I use the `rep_biclustermd` function for
this.

``` r
registerDoParallel(numCores) #'register' the number of CPUs for the loop

  l <- 0
 output <- foreach::foreach(l = seq_len(nrow(idx)), #use idx rows as loop counter
                             .combine = cbind, .packages = c('biclustermd')) %dopar% {
                              i <-  as.integer(idx$col_t[l]) #i is time/col clusters
                              j <-  as.integer(idx$row_s[l]) #j is space/row clusters 
                                
  #biclustering using i time and j site clusters:                             
    bc.tmp <- rep_biclustermd(bcdat, nrep = 500, col_clusters = i, row_clusters = j) 

            #get lowest sse across 500 reps
                    final.sse <- as.data.frame(bc.tmp$best_bc$SSE) 
                    final.sse <- min(final.sse$SSE)
                    
                    #calculate k:
            nsites*log(final.sse/nsites)+2*(i*j)  
    }
```

lets see what the `output` looks like:

``` r
dim(output)
```

    ## [1]   1 204

``` r
head(output[1,])
```

    ## result.1 result.2 result.3 result.4 result.5 result.6 
    ## 4240.206 4181.432 4161.682 4150.996 4142.501 4137.276

So we can see R outputted a single row and 204 columns, with the
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
value of each loop being stored in separate columns. The `foreach` loop
has labelled columns in the numeric order of the row index of the `idx`
dataframe. Since we know the order of rows of the `idx` dataframe is
congruent with the order of the output of the `foreach` loop, to
visualize the output I can prep a simple matrix with rows and columns
labelled according to the `idx` site and time clusters for plotting.

``` r
idx$out.k <- t(output)

#get a matrix for plotting:
mat.k <- matrix(idx$out.k, nrow=12, ncol=17)

#name rows and columns by cluster value
rownames(mat.k) <- c(2:13)
colnames(mat.k) <- c(2:18)

#get colours
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25) #PiYG

#plot
plot(mat.k, key=list(side=3, cex.axis=0.75), palette('Dark2'), digits=2, text.cell=list(cex=0.5), fmt.cell='%.0f',fmt.key="%.0f", col = coul, xlab = 'no. temporal clusters', ylab = 'no. spatial clusters', main = '', axis.row = 2)
```

![](AIC_EUstrat_12_10_2021.png)

The surface appears to be doing what is expected - I would expect to see
a continuous surface without major fluctuations in
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
values with incremental increases in cluster values. I might also expect
the lowest
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
values to be at an intermediate number of clusters, between under and
overfitting, depending on whether the range of clusters investigated was
selected appropriately. It looks like the optimal combination is
somewhere around 13 spatial and 11 temporal clusters. Let’s check this
properly, then run the biclustering using the optimal
![r.c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r.c "r.c")
combination:

``` r
#visualize a handful of best models
head(idx[order(idx$out.k),])
```

    ##     col_t row_s    out.k
    ## 120    11    13 3468.024
    ## 108    10    13 3471.626
    ## 96      9    13 3480.199
    ## 156    14    13 3480.934
    ## 119    11    12 3484.018
    ## 143    13    12 3484.812

Indeed, the best combination is 13 site and 11 time clusters, by &gt;2
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k")
values relative to the next best model. According to the AIC rule of
thumb, this means that it is unequivocally the best fitting model. So
let’s run it:

``` r
#get min k value
mn <- idx[idx$out.k == min(idx$out.k),] 

#run best biclustering
bestbc <- biclustermd::rep_biclustermd(bcdat, nrep = 500, col_clusters = mn$col_t, row_clusters = mn$row_s)
```

’Let’s see what the output looks like:

``` r
attributes(bestbc)
```

    ## $names
    ## [1] "best_bc" "rep_sse" "runtime"

The list called `best_bc` contains all the biclustering information
about the dataset and which site and time clusters each row and column
was assigned to. Let’s synthesize this information into a single
dataframe:

``` r
bestbc_assign <- gather(bestbc$best_bc)
head(bestbc_assign)
```

    ##   row_name   col_name row_cluster col_cluster bicluster_no     value
    ## 1        1 tsumm_1000           1           1            1 0.4226947
    ## 2       30 tsumm_1000           1           1            1 0.6996891
    ## 3       31 tsumm_1000           1           1            1 0.4221311
    ## 4       70 tsumm_1000           1           1            1 0.8536321
    ## 5       71 tsumm_1000           1           1            1 0.7292292
    ## 6       72 tsumm_1000           1           1            1 0.4701808

Each row of this dataframe represents a single cell of the input
dataset, showing its row and column names in the input dataset, and the
site and time clusters that the row was assigned to. This allows me to
get a list of the entity ID’s and their site cluster assignments, while
discarding the repetition of this information for each time point:

``` r
# get entity IDs and associated cluster:
map_clusters <- bestbc_assign[bestbc_assign$col_name == bestbc_assign$col_name[1], c("row_name", "row_cluster")]
length(unique(map_clusters$row_name));length(unique(map_clusters$row_cluster))
```

    ## [1] 1126

    ## [1] 13

In agreement with the input data and the biclustering specification
there are 1126 sites which have been put into and 13 site clusters.
Let’s name the columns appropriately:

``` r
names(map_clusters) <- c('ID_ENTITY', 'Site_cluster') 
head(map_clusters)
```

    ##   ID_ENTITY Site_cluster
    ## 1         1            1
    ## 2        30            1
    ## 3        31            1
    ## 4        70            1
    ## 5        71            1
    ## 6        72            1

Before we can visualize these clusters on a map, we need to attach the
spatial coordinates of each entity to this data. I use the input
dataframe which came with this data:

``` r
#get entity data:
ent <- forest[,c("lon","lat", "ID_ENTITY")]
ent <- unique(ent)

#merge to cluster data:
map_clusters <- merge(map_clusters, ent, by = 'ID_ENTITY')
head(map_clusters)
```

    ##   ID_ENTITY Site_cluster  lon  lat
    ## 1         1            1 -0.5 31.5
    ## 2        10            4 -0.5 40.5
    ## 3       100            1 -3.5 29.5
    ## 4      1000            6 29.5 37.5
    ## 5      1001            6 29.5 38.5
    ## 6      1002            6 29.5 39.5

The stratification is ready to be visualized. Let’s get a colour scheme
for the clusters and assign this to the grid cells for plotting.

``` r
#get spatial points with cluster values:

clust.points <- map_clusters[,!colnames(map_clusters) %in% c('sitenum', 'ID_ENTITY')]

#colour scheme:
colvar <- c('black','orange','forestgreen','chocolate3','blue','yellow4','hotpink','deepskyblue3','lawngreen','firebrick', 'skyblue','purple','coral4')  

#get a cluster colour dataframe and merge to spatial data:
coldf <- data.frame(colvar = colvar, Site_cluster = c(1:13))
clust.points <- merge(clust.points, coldf, by = 'Site_cluster')
head(clust.points)
```

    ##   Site_cluster  lon  lat colvar
    ## 1            1 -0.5 31.5  black
    ## 2            1 -5.5 30.5  black
    ## 3            1 -3.5 29.5  black
    ## 4            1 -9.5 30.5  black
    ## 5            1 -4.5 30.5  black
    ## 6            1 -7.5 32.5  black

I now have a dataframe of spatial points, their cluster assignments, and
a colour code representing the clusters. Let’s plot this:

``` r
#convert cluster values from numeric to factor, for plotting
clust.points$Site_cluster.ch <- as.character(clust.points$Site_cluster)
clust.points$Site_cluster.f <- factor(clust.points$Site_cluster.ch, levels = as.character(c(1:13)))

#plot:
ggplot() +
  geom_tile(data = clust.points, aes(x = lon, y = lat, fill = Site_cluster.f))+
  scale_fill_manual(values = unique(clust.points$colvar))+
  scale_x_continuous(limits = c(-11,44), breaks = seq(-10, 40, by=10)) +
  theme(legend.text=element_text(size=22),legend.position = 'bottom',legend.spacing.x = unit(0.5,'cm'),legend.spacing.y = unit(0.5,'cm'),legend.title=element_blank() , panel.background = element_rect(fill = "aliceblue"), axis.text = element_text(size=18, face='bold', colour = 'black'), axis.title = element_text(size=18, face='bold', colour = 'black'), title = element_text(size = 10, colour = 'black'))+
  scale_y_continuous(limits = c(29, 73), breaks = seq(30, 90, by=10)) +
  labs(x = "", y =  "")+
  geom_sf(data = world, fill="transparent", color="grey80")
```

![](biclustering_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

This is the final stratification representing spatial clusters of
coherent fire-relevant environmental change over the Holocene (numbered
and colour-coded in the legend). Some of the environmental patterns
agree with known environmental gradients, for example, the
high-elevation Swiss Alps and Moroccan Atlas Mountains are associated
with Scandinavia in cluster 9, highlighting the cold climates that these
regions share. This pattern is also seen in a modern European
environmental stratification generated by Metzger et al. (2005), as is
the spatial complexity of the clusters representing the Mediterranean
region, among other features. The similarity between the two
stratifications suggests the biclustering is successful at capturing
general environmental patterns native to Europe that have persisted over
the past 12000 years.

## References

J. Li, J. Reisner, H. Pham, S. Olafsson, and S. Vardeman. Biclustering
with missing data. *Information sciences*, 510:304–316, 2020.

Marc Joris Metzger, Robert Gerald Henry Bunce, Rob HG Jongman, Caspar A
Mucher, and John W Watkins. A climatic stratification of the environment
of europe. *Global* *ecology and biogeography*, 14(6):549–563, 2005.

A Mauri, BAS Davis, PM Collins, and Jed O Kaplan. The climate of europe
during the holocene: a gridded pollen-based reconstruction and its
multi-proxy evaluation. *Quaternary Science Reviews*, 112:109–127, 2015.

Marco Zanon, Basil AS Davis, Laurent Marquer, Simon Brewer, and Jed O
Kaplan. European forest cover during the past 12,000 years: a
palynological reconstruction based on modern analogs and remote sensing.
*Frontiers in plant science*, 9:253, 2018.
