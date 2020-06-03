from datetime import datetime, timedelta

from rise_set.angle import Angle
from rise_set.rates import ProperMotion
from rise_set.astrometry import calculate_airmass_at_times, make_ra_dec_target
from rise_set.visibility import Visibility

'''
    This example shows how to calculate the airmass values for a target at intervals throughout its visibility window
'''

def date_range_from_interval(start_time, end_time, dt=timedelta(minutes=15)):
    '''
        simple generator to generate datetimes evenly spaced throughout interval including start
    '''
    time = start_time
    while time < end_time:
        yield time
        time += dt


HOURS_PER_DEGREES = 15.0

# Maui telescope location details
site = {
    'latitude': Angle(degrees=20.7069444444),
    'longitude': Angle(degrees=-156.258055556),
    'horizon': Angle(degrees=25.0),
    'altitude': 3065,
    'ha_limit_neg': Angle(degrees=-4.533333 * HOURS_PER_DEGREES),
    'ha_limit_pos': Angle(degrees=4.4666667 * HOURS_PER_DEGREES)
}

# m23 target details
target = make_ra_dec_target(
    ra=Angle(degrees=269.267),
    dec=Angle(degrees=-18.985),
    ra_proper_motion=ProperMotion(Angle(degrees=0.000000328, units='arc')),
    dec_proper_motion=ProperMotion(Angle(degrees=-0.000000386, units='arc')),
    parallax=1.354,
    epoch=2000
)

# Extra data needed for computing visibility
airmass_limit = 1.6
moon_distance_limit = Angle(degrees=30.0)
start_date = datetime(2020, 6, 1)
end_date = datetime(2020, 6, 3, 23, 22, 0)

# Create a visibility object with the site parameters and a date range
visibility = Visibility(site, start_date, end_date, site['horizon'].in_degrees(), twilight='nautical', 
                        ha_limit_neg=-4.533333, ha_limit_pos=4.4666667, zenith_blind_spot=0)

# Gets the observable intervals for the target by intersecting the target_intervals, ha_intervals, and moon_distance intervals
observable_intervals = visibility.get_observable_intervals(target, airmass=airmass_limit, moon_distance=moon_distance_limit)

# chunk up the observable intervals into a list of time chunks (10 minutes) you want to compute the airmass at
chunked_intervals = []
for interval in observable_intervals:
    chunked_intervals.extend(
        [time for time in date_range_from_interval(interval[0], interval[1], dt=timedelta(minutes=10))])

# Compute the airmass at each datetime in the chunked up intervals
airmass_values = calculate_airmass_at_times(chunked_intervals, target, site['latitude'], site['longitude'], site['altitude'])

print("DATE - AIRMASS_VALUE")
for i, value in enumerate(chunked_intervals):
    print(f"{value} - {airmass_values[i]}")