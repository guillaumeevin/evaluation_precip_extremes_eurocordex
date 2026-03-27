import regionmask
import geopandas as gpd

# 1. Load the PRUDENCE regions
prudence = regionmask.defined_regions.prudence

# 2. Convert to GeoDataFrame
gdf = prudence.to_geodataframe()

# 3. Export to GeoJSON
gdf.to_file("regions.geojson", driver="GeoJSON")