# Cox-regressions
This is a script with cox regressions using the survival package
Known bug: test for proportionallity handles missing values in a weird way
If the last metabolite has no missing vaules, the loops seems to be working
This script has no normalization algorithm, instead, the values in the input are pre-normalized.
The script is created with two models for different adjustments
