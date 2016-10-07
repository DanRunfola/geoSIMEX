# geoSIMEX
geoSIMEX R Package

# Example of geoSIMEX
In R, open "Example/geoSIMEX_example.R"

# geoSIMEX documentation
In R, open "man/geoSIMEX.Rd" and click "preview"

# Data Requirements
geoSIMEX requires two datasets:
1. aid project data
2. country-level data (where the unit of analysis could be, for example, sub-counties). 

Both the aid finance and country-level data need to have hierarchies of administrative divisions. For example, in the aid finance data, if aid was allocated to a sub-county, information is needed about which district, region, and country that sub-county falls within. In country-level data, variables are needed for the sub-county, district, region and country. Names of these variables must be the same between the two datasets (e.g., if "ADM_3_NAME" is the variable name for sub counties in the country-level dataset, there must be a similar "ADM_3_NAME" in the aid project dataset.

# Needed Fixes
1. geoSIMEX needs ADM names to be entered as numbers; need to allow ADM info to be entered as strings
2. Something else



