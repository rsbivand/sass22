

## ---- echo=FALSE, eval=FALSE, results='hide'------------------------------------------------------------
## Population0 <-read.csv("data/KS101EW_oa11.csv")
## Population <- Population0[, 1:2]
## Employment0 <-read.csv("data/KS601EW_oa11.csv")
## Employment <- Employment0[, c(1:2, 6, 20)]
## all.equal(100*Employment[,3]/Employment[,2], Employment[,4])
## Employment <- merge(Employment, Population, by="GeographyCode")
## names(Employment) <- c("OA11CD", "all_categories_economic_activity", "economically_active_unemployed", "Unemployment", "Population")
## Sys.setenv(PROJ_NETWORK="ON")
## suppressPackageStartupMessages(library(sf))
## output_areas <- st_read("data/Camden_oa11.shp", quiet=TRUE)
## st_crs(output_areas) <- 27700
## oa_census <- merge(output_areas, Employment, by="OA11CD")
## st_write(oa_census, "oa_census.gpkg", delete_layer=TRUE, quiet=TRUE)


## -------------------------------------------------------------------------------------------------------
Sys.setenv(PROJ_NETWORK="ON")
library(sf)
packageVersion("sf")


## -------------------------------------------------------------------------------------------------------
camden <- st_read("oa_census.gpkg")


## -------------------------------------------------------------------------------------------------------
summary(with(camden, all_categories_economic_activity/Population))


## ---- message=FALSE-------------------------------------------------------------------------------------
ha <- st_area(camden)
library(units)
units(ha) <- as_units("ha")
dens <- camden[["all_categories_economic_activity"]]/ha
print(quantile(dens, seq(0, 1, 0.25)), digits=4)


## ---- results='hide', message=FALSE, warning=FALSE, cache=TRUE------------------------------------------
library(cartogram)
sf_cont <- cartogram_cont(camden, "all_categories_economic_activity", 7)


## -------------------------------------------------------------------------------------------------------
(cI <- classInt::classIntervals(camden[["Unemployment"]], n=6, style="fisher"))


## -------------------------------------------------------------------------------------------------------
pal <- rcartocolor::carto_pal((length(cI$brks)-1L), "TealGrn")


## ---- echo=TRUE, fig.width=7, fig.height=4, fig.cap="Figure 2: left panel: Percentage unemployment by output area in Camden, 2011 Census; right panel: cartogram of percentage unemployment with output areas adjusted to the size of the economically active population", eval=TRUE, fig=TRUE, message=FALSE----
library(tmap)
a <- tm_shape(camden) + tm_fill("Unemployment", breaks=cI$brks, palette=pal) + 
  tm_layout(main.title="Output Areas")
b <- tm_shape(sf_cont) + tm_fill("Unemployment", breaks=cI$brks, palette=pal) + 
  tm_layout(main.title="Cartogram", legend.show=FALSE)
tmap_arrange(a, b)


## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE-----------------------------------------------
## library(mapview)
## if (sf:::CPL_gdal_version() >= "3.1.0") mapviewOptions(fgb = FALSE)
## camden_wmap <- mapview(camden, zcol="Unemployment", col.regions=pal, at=cI$brks)
## mapshot(camden_wmap, file=file.path(getwd(), "camden_wmap.png"))


## -------------------------------------------------------------------------------------------------------
(degraded_proj4 <- st_crs(paste("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717",
"+x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs"))$proj4string)


## ---- echo=TRUE, eval=FALSE-----------------------------------------------------------------------------
## st_crs(camden)


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------
sfc <- st_centroid(camden[1,])


## -------------------------------------------------------------------------------------------------------
Sys.getenv("PROJ_NETWORK")


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------
sf_proj_network()


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------
(sfc_ll_best <- st_transform(sfc, "OGC:CRS84")) |> st_coordinates() |> print(digits=17)


## -------------------------------------------------------------------------------------------------------
trans_pipes <- sf_proj_pipelines(st_crs(sfc), "OGC:CRS84")
trans_pipes$accuracy


## -------------------------------------------------------------------------------------------------------
gsub("\\+step", "\n+step", trans_pipes$definition[10]) |> cat("\n")


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------
st_transform(sfc, "OGC:CRS84", pipeline=trans_pipes$definition[10]) |> st_distance(sfc_ll_best) |> c()


## -------------------------------------------------------------------------------------------------------
gsub("\\+step", "\n+step", trans_pipes$definition[11]) |> cat("\n")


## -------------------------------------------------------------------------------------------------------
st_transform(sfc, "OGC:CRS84", pipeline=trans_pipes$definition[11]) |> st_distance(sfc_ll_best) |> c()


## -------------------------------------------------------------------------------------------------------
gsub("\\+step", "\n+step", trans_pipes$definition[2]) |> cat("\n")


## -------------------------------------------------------------------------------------------------------
st_transform(sfc, "OGC:CRS84", pipeline=trans_pipes$definition[2]) |> st_distance(sfc_ll_best) |> c()


## ---- message=FALSE-------------------------------------------------------------------------------------
library(spdep)
(nb_q <- poly2nb(camden, queen=TRUE))
table(card(nb_q))


## -------------------------------------------------------------------------------------------------------
lw <- nb2listw(nb_q, style="W")


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------
sz1 <- object.size(lw)
sz2 <- object.size(listw2mat(lw))
requireNamespace("spatialreg", quietly=TRUE)
sz3 <- object.size(as(lw, "CsparseMatrix"))
sz <- c("listw"=sz1, "dense"=sz2, "sparse"=sz3)
round(sz/1024, 1)


## -------------------------------------------------------------------------------------------------------
moran.test(camden[["Unemployment"]], lw)


## -------------------------------------------------------------------------------------------------------
I_i <- localmoran(camden[["Unemployment"]], lw)
print(sum(I_i[,1])/Szero(lw), digits=9)


## -------------------------------------------------------------------------------------------------------
I_ia <- localmoran(camden[["Unemployment"]], lw, mlvar=FALSE)
print(sum(I_ia[,1])/Szero(lw), digits=16)


## ---- echo=TRUE, fig.width=7, fig.height=4, fig.cap="Figure 5: Local Moran's I hotspot cluster cores: left panel: choropleth map; right panel: cartogram", eval=TRUE, fig=TRUE----
HS <- attr(I_ia, "quadr")$mean
is.na(HS) <- p.adjust(I_ia[, "Pr(z != E(Ii))"], "none") >= 0.05
camden$HS <- HS
sf_cont$HS <- HS
a <- tm_shape(camden) + tm_fill("HS", palette=pal, colorNA="grey95", textNA="Not \"interesting\"") + tm_layout(main.title="Output Areas")
b <- tm_shape(sf_cont) + tm_fill("HS", palette=pal, colorNA="grey95", textNA="Not \"interesting\"") + tm_layout(main.title="Cartogram", legend.show=FALSE)
tmap_arrange(a, b)

