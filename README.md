# euclid_phz_testing

This repositiory contains a set of tools to run and analyse Euclid photo-z data challenges.

## photoz_metrics.py

Script to compute the photo-z performance metrics. The `scatter` option plots the photo-z value against a reference redshift and computes the usual photo-z statistics (such as sigma, outlier rate and bias). The `PDF` option plots the n(z), PDF(z-z_ref) and statistics.

```
usage: photoz_metrics.py [-h] [-i INPUT] [-o OUTPUT] [-input_type INPUT_TYPE]
                         [-zmin ZMIN] [-zmax ZMAX] [-title TITLE]
                         [-select SELECT] [-density]
                         [-stats_output STATS_OUTPUT] [-z_bins Z_BINS]
                         task

positional arguments:
  task                  Metric to compute and plot. [scatter, nz_bins, PDF,
                        bias, all]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -o OUTPUT, --output OUTPUT
                        output file
  -input_type INPUT_TYPE
                        input file type (default: ECLD_PHZ)
  -zmin ZMIN            minimum redshift for the scatter plot (default: 0.0)
  -zmax ZMAX            maximum redshift for the scatter plot (default: 6.0)
  -title TITLE          title (default: None)
  -select SELECT        selection string (default: None)
  -density              display scatter plot as density
  -stats_output STATS_OUTPUT
                        stats output file
  -z_bins Z_BINS        redshift bins for PDF analysis (default:
                        0.2,0.45,0.55,0.7,0.8,0.9,1.0,1.15,1.35,1.65,2.0)
```

## simulate.py

Script to simulate fluxes and flux errors, given AB magnitude
depths and sky backgrounds.

```
usage: simulate.py [-h] [-i INPUT] [-o OUTPUT] [-seed SEED] [-r_ref R_REF]
                   [-filter_names FILTER_NAMES]
                   [-filter_skies_AB FILTER_SKIES_AB]
                   [-filter_depths_AB FILTER_DEPTHS_AB] [-repeat REPEAT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -o OUTPUT, --output OUTPUT
                        output file
  -seed SEED            random seed
  -r_ref R_REF          reference source size on the sky in which depths are defined
  -filter_names FILTER_NAMES
                        filter names
  -filter_skies_AB FILTER_SKIES_AB
                        sky brightness (AB) in filters
  -filter_depths_AB FILTER_DEPTHS_AB
                        filter depths (AB)
  -repeat REPEAT        Number of repetitions
  ```

## simulate_from_inst.py

Script to simulate fluxes and flux errors, given AB magnitude
depths and sky backgrounds, from true instrumental
characteristics.


```
usage: simu_obs.py [-h] [-i INPUT] [-o OUTPUT] [-seed SEED] option

positional arguments:
  option                Which action

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file(s)
  -o OUTPUT, --output OUTPUT
                        output file(s)
  -seed SEED            Random seed. Default: 20091982
```

## get_dust_corr.py

script to load dust map and get extinction value and/or apply corrections.

Some typical coefficients useful for Euclid:

* MegaCam-uS: 4.018
* MegaCam-u: 4.138
* HSC-G: 3.258
* HSC-R: 2.286
* HSC-I: 1.641
* HSC-Z: 1.264
* HSC-Y: 1.081
* DECam-G: 3.687
* DECam-R: 2.499
* DECam-I: 1.873
* VIS: 1.780
* DECam-Z: 1.399
* VISTA-Y: 1.167
* VISTA-J: 0.965
* VISTA-H: 0.637


```
usage: get_dust_corr.py [-h] [-c COLUMNS] [-s SOUTHFILE] [-n NORTHFILE]
                        [-band BAND] [-coef COEF] [-correct] [-add_corr]
                        input output

positional arguments:
  input                 input file
  output                output file

optional arguments:
  -h, --help            show this help message and exit
  -c COLUMNS, --columns COLUMNS
                        Names of coordinates colums (default: ra,dec)
  -s SOUTHFILE, --southFile SOUTHFILE
                        Dust map (south galactic)
  -n NORTHFILE, --northFile NORTHFILE
                        Dust map (north galactic)
  -band BAND            List of bands to correct the extinction for
  -coef COEF            Corresponding Albda/E(B-V) (in mags)
  -correct              Correct fluxes for extinction
  -add_corr             Add flux correction factor for extinction
```
