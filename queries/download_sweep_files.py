# download_sweep_files.py
import math
import urllib


# Use the eFEDS field as an example
ra_min, ra_max = 126.5, 145.5
dec_min, dec_max = -3.2, 6.2


# Right ascension is taken in steps of 10, declination in steps of five
ra_min = 10 * math.floor(ra_min / 10)
ra_max = 10 * math.ceil(ra_max / 10)
dec_min = 5 * math.floor(dec_min / 5)
dec_max = 5 * math.ceil(dec_max / 5)


def sgnstr(dec):
    """Returns 'p' if dec is positive, 'm' if it's negative."""
    return "p" if dec >= 0 else "m"


def give_region_string(ra, dec):
    """Returns a string conforming to the naming convention of the SWEEP
    catalogues following the <AAA>c<BBB> pattern where AAA is the ra, c 
    the p or m for the sign of dec and BBB is the dec."""
    return f"{ra:03}{sgnstr(dec)}{abs(dec):03}"


web_base = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/"
count = (ra_max - ra_min) / 10 * abs(dec_max - dec_min) / 5
print(f"Fetching {int(count)} regions with ra between {ra_min} and {ra_max},"
      f" dec between {dec_min} and {dec_max}")
directory = "sweep/"  # Directory to store the data in
i = 1
for ra in range(ra_min, ra_max, 10):
    for dec in range(dec_min, dec_max, 5):
        reg_min, reg_max = give_region_string(ra, dec), \
            give_region_string(ra + 10, dec + 5)
        region = f"sweep-{reg_min}-{reg_max}.fits"
        website = web_base + region
        try:
            site = urllib.request.urlopen(website)
            print(f"Downloading {region}... Filesize:"
                  f" {round((site.length/1024**3),2)} GB")
            urllib.request.urlretrieve(website, f"{directory}/SWEEP/{region}")
            print(
                f"Downloaded and stored {region} in {directory} ({i} of {int(count)})")
        except urllib.error.HTTPError:
            print(f"Couldn't find {region}. Moving on to the next file.")
        i += 1
