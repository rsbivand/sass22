import geopandas
geopandas.show_versions()


# avoid https://github.com/Toblerity/Fiona/issues/986
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

geodf = geopandas.read_file('oa_census.gpkg')
cent1 = geodf.centroid[0]
cx = cent1.x
cy = cent1.y
geodf.crs

import pyproj
from pyproj.transformer import TransformerGroup
pyproj.network.set_network_enabled(True)
tg = TransformerGroup(geodf.crs, 'OGC:CRS84', always_xy=True)
itg = TransformerGroup('OGC:CRS84', geodf.crs, always_xy=True)
tg

import re
print(re.sub("step", "\nstep", tg.transformers[0].definition))

ll0 = tg.transformers[0].transform(cx, cy)
ll0

print(re.sub("step", "\nstep", tg.transformers[2].definition))
ll2 = tg.transformers[2].transform(cx, cy)

from shapely.geometry import Point
point_df = geopandas.GeoDataFrame({'geometry': [Point(itg.transformers[0].transform(ll2[0], ll2[1]))]}, crs='EPSG:27700')
point_df.distance(cent1)

print(re.sub("step", "\nstep", tg.transformers[10].definition))
ll10 = tg.transformers[10].transform(cx, cy)
point_df = geopandas.GeoDataFrame({'geometry': [Point(itg.transformers[0].transform(ll10[0], ll10[1]))]}, crs='EPSG:27700')
point_df.distance(cent1)

import numpy as np
dens = geodf['all_categories_economic_activity'] / (geodf.area/10000)
arr = np.quantile(dens, np.linspace(0, 1, 5))
np.set_printoptions(precision=3, suppress=True)
print(arr)

from cartogram_geopandas import make_cartogram
transformed_geodf = make_cartogram(geodf, 'all_categories_economic_activity', 7, inplace=False)

import mapclassify as mc
fj6 = mc.FisherJenks(geodf["Unemployment"], k=6)
fj6

cIfj = mc.UserDefined(geodf["Unemployment"], bins=[2.027740,  3.668135,  5.542755,  7.745236, 11.191642, 18.623482])
cIfj

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
mycmap = ListedColormap(["#b0f2bc", "#82e6aa", "#5bd4a4", "#3fbba3", "#2e9ea1", "#257d98"])
plt.figure(figsize=(13,7))
ax1 = plt.subplot(121)
cIfj.plot(geodf, cmap=mycmap, ax=ax1, legend=True, legend_kwds={'loc':'lower left'});
plt.title('Output Areas')
plt.axis('off');
ax2 = plt.subplot(122)
cIfj.plot(transformed_geodf, cmap=mycmap, ax=ax2);
plt.title('Cartogram')
plt.axis('off');
plt.show()

from libpysal import weights
w_cont = weights.Queen.from_dataframe(geodf)
w_cont.histogram

import esda
w_cont.transform = "R"
mi = esda.Moran(geodf["Unemployment"], w_cont)
np.set_printoptions(precision=9, suppress=True)
print(np.array((mi.I, mi.z_rand)))

lisa = esda.Moran_Local(geodf["Unemployment"], w_cont)
np.sum(lisa.Is)/w_cont.s0

