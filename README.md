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
  task                  Metric to compute and plot [scatter, PDF]

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

Script to simulate the flux and flux errors, given AB magnitude depths and sky backgrounds.

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
