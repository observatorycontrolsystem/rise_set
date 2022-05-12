from datetime import datetime, timedelta

from rise_set.angle import Angle
from rise_set.rates import ProperMotion
from rise_set.visibility import Visibility
from rise_set.astrometry import make_ra_dec_target

"""
    This example shows how to calculate the observable intervals of a target for a given site.
    It also shows how to calculate the intermediate visibility windows if they are useful.
"""

HOURS_PER_DEGREES = 15.0

# Maui telescope location details
site = {
    'latitude': Angle(degrees=20.7069444444),
    'longitude': Angle(degrees=-156.258055556),
    'horizon': Angle(degrees=25.0),
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
maximum_moon_phase = 1.0
start_date = datetime(2020, 6, 1)
end_date = datetime(2020, 6, 3, 23, 22, 0)

# Create a visibility object with the site parameters and a date range
visibility = Visibility(site, start_date, end_date, site['horizon'].in_degrees(), twilight='nautical', 
                        ha_limit_neg=-4.533333, ha_limit_pos=4.4666667, zenith_blind_spot=0)

# Get Dark intervals to get the sun up / sun down intervals for a given site
dark_intervals = visibility.get_dark_intervals()
print(f"dark_intervals = {dark_intervals}")

# Get Target Intervals for set of intervals when target is above horizon for the site
target_intervals = visibility.get_target_intervals(target, up=True, airmass=airmass_limit)
print(f"target_intervals = {target_intervals}")

# Get the HA Intervals for the set of intervals when the target is within the telescopes hour angle limits
ha_intervals = visibility.get_ha_intervals(target)
print(f"ha_intervals = {ha_intervals}")

# Get the moon distance intervals over which the target is above the horizon and not within moon_distance_limit of the moon
moon_distance_intervals = visibility.get_moon_distance_intervals(target, target_intervals, moon_distance=moon_distance_limit,
                                                                 chunk_size=timedelta(minutes=30))
print(f"moon_distance_intervals = {moon_distance_intervals}")

# Get the moon phase intervals over which either the moon is below the horizon or the phase is less than maximum_moon_phase
moon_phase_intervals = visibility.get_moon_phase_intervals(target_intervals, max_moon_phase=maximum_moon_phase,
                                                           chunk_size=timedelta(minutes=30))
print(f"moon_phase_intervals = {moon_phase_intervals}")

# Gets the observable intervals for the target by intersecting the target_intervals, ha_intervals, and moon_distance, and moon phase intervals
observable_intervals = visibility.get_observable_intervals(target, airmass=airmass_limit, moon_distance=moon_distance_limit, moon_phase=maximum_moon_phase)
print(f"observable_intervals = {observable_intervals}")
