# weave-streams
Target selection for stellar streams in WEAVE 

# Usage

```
> make_widestreams_catalogue [gd1, pal5, orphan, tripsc]
> make_targetlist [gd1, pal5, orphan, tripsc] --input <dataframes produced in previous step>
```

The first step will preselect the data in position from the full Gaia eDR3
catalogue xmatched with Pan-STARRS DR1. It will output a dataframe (hdf5) with
position selected catalogue plus boolean columns corresponding to the different
filters it computed (color-magnitude, proper motion, quality, etc).

The second step will take these dataframes and apply the filters, save the
filtered dataframes, and concatenate all in a single FITS file. Any external catalogue
[sag, cetus] will also be added to this.

To run this you need the full Gaia eDR3 as a single hdf5 on disk and `vaex` to
handle it outside of memory. Alternatively you can write your own data query
routine and get the data from elsewhere.

One general configuration file is necessary (default is `config.yaml`; see
examples). All other configuration is internal and can be found in the `conf`
directory. 

## Get set of isocrhones 

get it from Zenodo a set of MIST 2.1 isochrones using:
```
curl --cookie-jar zenodo-cookies.txt "https://zenodo.org/record/6600715/files/isocs.tar.gz?download=1" | tar -xvz
```
