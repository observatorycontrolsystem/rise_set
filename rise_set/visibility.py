"""``visibility.py`` - Visibility class for calculating the observable intervals.

Examples:
    Get the observable intervals of a target:

    >>> from rise_set.angle import Angle
    >>> from rise_set.visibility import Visibility
    >>> visibility = Visibility(
    >>>     site_details, start_date, end_date, horizon, twilight='nautical'
    >>> )
    >>> observable_intervals = visibility.get_observable_intervals(
    >>>     target, airmass=1.6, moon_distance=Angle(degrees=30), moon_phase=0.9
    >>> )
"""

# Required for true (non-integer) division
from __future__ import division
from builtins import range
from builtins import object

# Standard libary imports
import datetime
import math
import copy

# Internal imports
from rise_set.astrometry import (
    calc_local_hour_angle, calc_sunrise_set, calc_planet_rise_set, calc_rise_set, RiseSetError,
    Star, gregorian_to_ut_mjd, ut_mjd_to_gmst, date_to_tdb, apparent_planet_pos, calculate_moon_phase,
    calculate_zenith_distance, mean_to_apparent, angular_distance_between, elem_to_topocentric_apparent)
from rise_set.angle import Angle
from rise_set.exceptions import InvalidHourAngleLimit
from rise_set.moving_objects import find_moving_object_up_intervals
from rise_set.utils import (
    coalesce_adjacent_intervals, intersect_intervals, is_sidereal_target, inverse_intervals,
    intersect_many_intervals, is_moving_object, is_static_target, target_to_jform)

# Import logging modules
import logging

# Set the logger name
_log = logging.getLogger('rise_set.visibility')


# Set convenient constants
ONE_DAY  = datetime.timedelta(days=1)
MIDNIGHT = datetime.time()
# The average moon refraction from Astronomical Almanac A12: 34 minutes
MOON_REFRACTION = Angle(degrees=-0.5666667)


def set_airmass_limit(airmass, horizon):
    ''' Compare the provided maximum airmass limit with the horizon of a telescope, and
        return the effective horizon for rise/set purposes (whichever is higher elevation). If
        no airmass is provided, we default to the horizon.

        airmass = 1 / cos(zenith)
        horizon = 90 - zenith'''

    if not airmass:
        return horizon

    # We convert the horizon to airmass, not vice versa, to avoid small angle problems if a
    # very large airmass is provided
    zenith_distance = 90 - horizon
    horizon_airmass = 1 / math.cos(math.radians(zenith_distance))

    effective_horizon = horizon
    if airmass < horizon_airmass:
        effective_horizon = 90 - math.degrees(math.acos(1 / airmass))

    return effective_horizon


class Visibility(object):

    """The Visibility class is used to calculate target visibilities for a given site.

    The Visibility class is instantiated with a given site and time range. 
    It can be used repeatedly with different targets to get those targets 
    observable intervals at the given site. It takes ha_limits, horizon limits, 
    airmass limits, moon angular distance limits, and zenith blind spot limits 
    into account when computing over all observable intervals for a target.

    Args:
        site (dict): Dictionary of site properties. Should contain Angles for latitude and longitude
        start_date (datetime): Start datetime over which you want to calculate intervals
        end_date (datetime): End datetime over which you want to calculate intervals
        horizon (float): Horizon angle in degrees for the site
        twilight (str): Type of twilight to use for rise/set calculation. Can be one of `sunrise`, `sunset`, `civil`, `nautical`, `astronomical`
        ha_limit_neg (float): The hour angle negative limit for the telescope
        ha_limit_pos (float): The hour angle positive limit ror the telescope
        zenith_blind_spot (float): blind spot angle in degrees over which the telescope cannot observe
    Raises:
        rise_set.exceptions.InvalidHourAngleLimit: If the positive or negative hour angles provided are out of the possible range.
    """
    def __init__(self, site, start_date, end_date, horizon=0, twilight='sunrise',
                 ha_limit_neg=-4.9, ha_limit_pos=4.9, zenith_blind_spot=0):
        self.site         = site
        self.start_date   = start_date
        self.end_date     = end_date
        self.horizon      = Angle(degrees=horizon)
        self.twilight     = twilight
        self.zenith_blind_spot = Angle(degrees=zenith_blind_spot)

        if ha_limit_pos > 12.0 or ha_limit_pos < 0.0:
            msg = "Positive hour angle limit must fall between 0 and 12 hours"
            raise InvalidHourAngleLimit(msg)

        if ha_limit_neg < -12.0 or ha_limit_neg > 0.0:
            msg = "Negative hour angle limit must fall between -12 and 0 hours"
            raise InvalidHourAngleLimit(msg)

        self.ha_limit_neg = ha_limit_neg
        self.ha_limit_pos = ha_limit_pos

        self.dark_intervals = []
        self.moon_dark_intervals = []


    def get_dark_intervals(self):
        """Returns the dark intervals for the site.
        
        Returns the night time dark intervals for the given site and date range set
        in this visibility object.

        Returns:
            list: A list of tuples of start/end datetime pairs that make up the dark intervals for this site.
        """
        # Don't compute this again if we've already done it
        if self.dark_intervals:
            return self.dark_intervals

        target = 'sun'

        self.dark_intervals = self.get_target_intervals(target, up=False)

        return self.dark_intervals


    def get_moon_dark_intervals(self):
        """Returns the dark moon intervals for the site.
        
        Returns the dark moon intervals (moon is not visible) for the given
        site and date range set in this visibility object.

        Returns:
            list: A list of tuples of start/end datetime pairs that make up the dark moon intervals for this site.
        """
        # Don't compute this again if we've already done it
        if self.moon_dark_intervals:
            return self.moon_dark_intervals

        target = 'moon'

        self.moon_dark_intervals = self.get_target_intervals(target, up=False)

        return self.moon_dark_intervals

    def _add_moon_interval(self, time, target_app_ra, target_app_dec, constraint=Angle(degrees=30)):
        moon_app_ra, moon_app_dec, diameter = apparent_planet_pos('moon', time['tdb'], self.site)
        # call slalib to get the angular moon distance
        target_moon_dist = angular_distance_between(target_app_ra, target_app_dec, moon_app_ra, moon_app_dec)
        # if that moon distance is > the constraint, add this interval to final intervals
        return target_moon_dist.in_degrees() >= constraint.in_degrees()

    def _add_zenith_interval(self, time, target_app_ra, target_app_dec, constraint=Angle(degrees=0)):
        ha = calc_local_hour_angle(target_app_ra, self.site['longitude'], time['time'])
        target_zenith_dist = calculate_zenith_distance(self.site['latitude'], target_app_dec, ha)
        # if that zenith distance is > the constraint, add this interval to final intervals
        return target_zenith_dist.in_degrees() >= constraint.in_degrees()

    def _get_chunked_intervals(self, target, target_intervals, compare_func, constraint, chunk_size=datetime.timedelta(minutes=30)):
        '''Returns a set of datetime 2-tuples, each of which represents an interval
           of time that the target is greater than the constraint away from the thing-to-be-avoided.
           The supplied compare_func calculates the distance to it's specific obstacle (moon, zenith).
        '''
        intervals = []

        for start, end in target_intervals:
            chunk_start = start
            chunk_end = min(chunk_start + chunk_size, end)
            while chunk_start != chunk_end and chunk_end <= end:
                # get the tdb date of the start time of the interval
                tdb = date_to_tdb(chunk_start)
                # get the apparent ra/dec for the target, and for the moon at this timestamp
                if is_sidereal_target(target):
                    target_app_ra, target_app_dec = mean_to_apparent(target, tdb)
                else:
                    target_app_ra, target_app_dec = elem_to_topocentric_apparent(chunk_start, target, self.site,
                                                                                 target_to_jform(target))

                if compare_func({'time': chunk_start, 'tdb': tdb}, target_app_ra, target_app_dec, constraint):
                    intervals.append((chunk_start, chunk_end))

                # increment the chunk_start/end up
                chunk_start = chunk_end
                chunk_end = min(chunk_start + chunk_size, end)

        intervals = coalesce_adjacent_intervals(intervals)
        return intervals

    def get_moon_distance_intervals(self, target, target_intervals, moon_distance=Angle(degrees=30), chunk_size=datetime.timedelta(minutes=30)):
        """Returns the moon distance intervals for the given target.
        
        Returns the intervals for which the given target is greater than the angular moon distance
        away from the moon at the given site and date range set in this visibility object.

        Args:
            target (dict): A dictionary of target details in the rise-set library format
            target_intervals (list): A list of datetime tuples that represent the above horizon intervals for the target. Returned by get_target_intervals()
            moon_distance (Angle): The minimum angular moon distance that the target must be away from the moon
            chunk_size (timedelta): The time delta over which to calculate if the target intervals are out of range of the zenith.
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is greater than angular moon distance away from moon.
        """
        return self._get_chunked_intervals(target, target_intervals, self._add_moon_interval, moon_distance, chunk_size)

    def get_zenith_distance_intervals(self, target, target_intervals, chunk_size=datetime.timedelta(minutes=1)):
        """Returns the zenith distance intervals for the given target.

        Returns the intervals for which the given target is greater than zenith distance away from the zenith at the given site 
        and date range set in this visibility object.

        Args:
            target (dict): A dictionary of target details in the rise-set library format
            target_intervals (list): A list of datetime tuples that represent the above horizon intervals for the target.
                                     Returned by get_target_intervals()
            chunk_size (timedelta): The time delta over which to calculate if the target intervals are out of range of
                                    the zenith.
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is greater than zenith distance away from zenith.
        """
        return self._get_chunked_intervals(target, target_intervals, self._add_zenith_interval, self.zenith_blind_spot, chunk_size)

    def get_moon_phase_intervals(self, target_intervals, max_moon_phase=1.0, chunk_size=datetime.timedelta(minutes=30)):
        """Returns the intervals in which the target is visible and the moon phase does not exceed moon_phase

        Args:
            target_intervals (list): A list of datetime tuples that represent the above horizon intervals for the target.
                                     Returned by get_target_intervals()
            moon_phase (float): max allowablable fractional moon phase from 0 to 1 (full moon).
            chunk_size (timedelta): The time delta over which to calculate if the target intervals are out of range of
                                   the zenith.
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is visible and the moon phase is less than max_moon_phase
        """
        # If target intervals are empty, no point in calculating moon phases
        if not target_intervals:
            return target_intervals
        bad_intervals = []
        moon_up_intervals = self.get_target_intervals(target="moon", up=True)
        target_moon_up_intervals = intersect_intervals(target_intervals, moon_up_intervals)
        # Here we collect the bad intervals where the moon is up and the moon phase is greater than the max constraint
        for start, end in target_moon_up_intervals:
            chunk_start = start
            chunk_end = min(chunk_start + chunk_size, end)
            while chunk_start != chunk_end and chunk_end <= end:
                moon_phase = calculate_moon_phase(chunk_start, self.site['latitude'].in_radians(), self.site['longitude'].in_radians())
                if moon_phase > max_moon_phase:
                    bad_intervals.append((chunk_start, chunk_end))

                # increment the chunk_start/end up
                chunk_start = chunk_end
                chunk_end = min(chunk_start + chunk_size, end)

        bad_intervals = coalesce_adjacent_intervals(bad_intervals)
        # Bad intervals are inversed to get good intervals and then intersected with target intervals to get target intervals where moon_phase constraint is met
        good_intervals = inverse_intervals(bad_intervals, target_intervals[0][0], target_intervals[-1][1])
        good_intervals = intersect_intervals(good_intervals, target_intervals)
        return good_intervals

    def get_target_intervals(self, target, up=True, airmass=None):
        """Returns the above or below horizon intervals for the given target.
        
        Returns the above (up=True) or below (up=False) horizon intervals for the given target and given site and date range set in this visibility object.

        Args:
            target (dict): A dictionary of target details in the rise-set library format
            airmass (float): The maximum acceptable airmass for this target to be observable in
            up (boolean): True (default) if you want intervals above the horizon, False for below
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is above/below the horizon.
        Raises:
            rise_set.exceptions.RiseSetError: If there was a problem calculating the rise/set/transfer times of the target
        """
        effective_horizon = set_airmass_limit(airmass, self.horizon.in_degrees())

        # Return dark intervals for static targets
        if is_static_target(target):
            intervals = self.get_dark_intervals()
        # Handle moving objects differently from stars
        elif is_moving_object(target):
            intervals = self._get_moving_object_target_intervals(target, effective_horizon)
        # The target has an RA/Dec
        else:
            intervals = self._get_ra_target_intervals(target, up, airmass, effective_horizon)

        return intervals

    def _get_moving_object_target_intervals(self, target, effective_horizon):
        window = {
                   'start' : self.start_date,
                   'end'   : self.end_date,
                 }
        site = self.site.copy()
        site['horizon'] = Angle(degrees=effective_horizon)
        intervals, _ = find_moving_object_up_intervals(window, target, site)

        return intervals


    def _get_ra_target_intervals(self, target, up, airmass, effective_horizon):
        star = Star(self.site['latitude'], target, effective_horizon)

        if up:
            day_interval_func = self._find_when_target_is_up
        else:
            day_interval_func = self._find_when_target_is_down

        # Find rise/set/transit for each day
        intervals = []
        current_date = self.start_date
        while current_date < self.end_date + ONE_DAY:
            one_day_intervals = day_interval_func(target, current_date, star, airmass)

            # Add today's intervals to the accumulating list of intervals
            intervals.extend(one_day_intervals)

            # Move on to tomorrow
            current_date += ONE_DAY

        # Collapse adjacent intervals into continuous larger intervals
        intervals = coalesce_adjacent_intervals(intervals)
        intervals = intersect_intervals(intervals, [(self.start_date,self.end_date)])

        return intervals


    def get_ha_intervals(self, target):
        """Returns the hour angle intervals for the given target.

        Returns the hour anle intervals for the given target and given site and date range set in this visibility object.
        The hour angle intervals are uninterupted chunks of time that the target is within the hour angle limits of the
        telescope.

        Args:
            target (dict): A dictionary of target details in the rise-set library format
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is within HA limits.
        """
        SIDEREAL_SOLAR_DAY_RATIO = 1.002737909350
        SIDEREAL_SOLAR_DAY = datetime.timedelta(seconds=(ONE_DAY.total_seconds() / SIDEREAL_SOLAR_DAY_RATIO))

        earliest_date = self.start_date - SIDEREAL_SOLAR_DAY

        tdb    = date_to_tdb(earliest_date)
        mjd    = gregorian_to_ut_mjd(earliest_date)
        gmst   = ut_mjd_to_gmst(mjd)

        # Need the apparent ra/dec for getting correct ha limits on high dec targets
        apparent_ra, apparent_dec = mean_to_apparent(target, tdb)

        # Flip the neg/pos ha limits if site is in the southern hemisphere
        ha_neg = self.ha_limit_neg
        ha_pos = self.ha_limit_pos
        if self.site['latitude'].in_degrees() < 0:
            ha_neg = -self.ha_limit_pos
            ha_pos = -self.ha_limit_neg

        # the rise time
        hour_rise = ha_neg + apparent_ra.in_hours() - \
            self.site['longitude'].in_hours() - gmst.in_hours()
        hour_rise /= SIDEREAL_SOLAR_DAY_RATIO

        # the set time
        hour_set  = ha_pos + apparent_ra.in_hours() - \
            self.site['longitude'].in_hours() - gmst.in_hours()
        hour_set /= SIDEREAL_SOLAR_DAY_RATIO

        current_rise = earliest_date + datetime.timedelta(hours=hour_rise)
        current_set = earliest_date + datetime.timedelta(hours=hour_set)

        # Find hour angle limits for each day
        intervals = []
        while current_set < (self.end_date + SIDEREAL_SOLAR_DAY):
            intervals.append((current_rise, current_set))
            current_rise += SIDEREAL_SOLAR_DAY
            current_set += SIDEREAL_SOLAR_DAY

        # do not exceed start/end dates
        intervals = coalesce_adjacent_intervals(intervals)
        intervals = intersect_intervals(intervals,[(self.start_date,self.end_date)])

        return intervals


    def get_observable_intervals(self, target, airmass=None, moon_distance=Angle(degrees=30), moon_phase=1.0):
        """Returns the observable intervals for the given target.
        
        Returns the observable intervals for the given target and given site and date range set in this visibility object.
        The observable intervals are the intersections of the dark intervals at the site, the above horizon target intervals,
        the moon distance intervals of the target, the hour angle intervals of the target, and the zenith blind spot intervals
        of the target.

        Args:
            target (dict): A dictionary of target details in the rise-set library format
            airmass (float): The maximum acceptable airmass for this target to be observable in
            moon_distance (Angle): The minimum acceptable angular distance between the moon and the target
            moon_phase (float): The maximum acceptable moon phase fraction from 0 to 1 (full moon)
        Returns:
            list: A list of tuples of start/end datetime pairs that make up the intervals over which this target is observable.
        Raises:
            rise_set.exceptions.RiseSetError: If there was a problem calculating the rise/set/transfer times of the target
        """
        # get the intervals of each separately
        dark               = self.get_dark_intervals()
        above_horizon      = self.get_target_intervals(target, airmass=airmass)

        # Apply the moon distance constraint
        if moon_distance.in_degrees() <= 0.5 or is_static_target(target):
            moon_avoidance = above_horizon
        else:
            moon_avoidance = self.get_moon_distance_intervals(target, above_horizon, moon_distance)

        # Apply the moon phase constraint
        if moon_phase >= 0.99:
            # If its within 1% of 1.0, ignore the constraint and just return all intervals
            moon_constraints = moon_avoidance
        else:
            moon_constraints = self.get_moon_phase_intervals(moon_avoidance, moon_phase)

        if is_sidereal_target(target):
            within_hour_angle = self.get_ha_intervals(target)
        else:
            # if the target type is such that there is no 'ra'/'dec' values, then we cannot calculate the ha intervals,
            # so we just use the target intervals instead (since they are all intersected together next). This is true
            # for moving objects and static objects.
            within_hour_angle = above_horizon

        if self.zenith_blind_spot.in_degrees() <= 0.0 or is_static_target(target):
            zenith_hole_avoidance = above_horizon
        else:
            zenith_hole_avoidance = self.get_zenith_distance_intervals(target, above_horizon)

        # find the overlapping intervals between them
        intervals = intersect_many_intervals(dark, above_horizon, within_hour_angle,
                                             moon_constraints, zenith_hole_avoidance)

        return intervals


    def _find_when_target_is_down(self, target, dt, star=None, airmass=None):
        '''Returns a single datetime 2-tuple, representing an interval
           of uninterrupted time below the horizon at the specified site, for the
           requested date.

           Note: Even though this function currently ignores times, the dt object
           must be a datetime, *not* a date.
        '''

        # Ensure we only deal with dates, because our rise/set/transit tuple is
        # day specific.
        # TODO: Extend to arbitrary start/end times
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

        # We will calculate down intervals as the inverse of the up intervals
        _log.debug("dt: %s", dt)
        up_intervals = self._find_when_target_is_up(target, dt, star, airmass)

        if not up_intervals:
            _log.warn("Got no up intervals!")
            _log.warn("dt was: %s", dt)
            _log.warn("target was: %s", target)

        down_intervals = []

        # If the first value has time 00:00:00, then the target starts up
        if up_intervals[0][0].time() == MIDNIGHT:
            pass

        # Otherwise the target starts down - so there's one extra interval at start
        else:
            down_start = dt
            down_end   = up_intervals[0][0]

            down_intervals.append((down_start, down_end))


        # Proceed through the intervals, extracting the gaps
        for i in range(len(up_intervals) - 1):
            down_start = up_intervals[i][1]
            down_end   = up_intervals[i+1][0]

            down_intervals.append((down_start, down_end))


        # If the target sets before the end of the day, grab that as an
        # extra down interval
        if up_intervals[-1][1].time() != MIDNIGHT:
            down_start = up_intervals[-1][1]
            down_end   = dt + ONE_DAY

            down_intervals.append((down_start, down_end))

        return down_intervals


    def _find_when_target_is_up(self, target, dt, star=None, airmass=None):
        '''Returns a single datetime 2-tuple, representing an interval
           of uninterrupted time above the horizon at the specified
           site, for the requested date.

           Note: Even though this function currently ignores times, the dt
           object must be a datetime, *not* a date.

           TODO: Clean up this unpleasant interface - no star need be passed if
           the sun is the target, but one *must* be passed otherwise. This isn't
           obvious from the method signature.
        '''
        effective_horizon = set_airmass_limit(airmass, self.horizon.in_degrees())

        # Remove any time component of the provided datetime object
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

        # Get the rise/set/transit times for the target, for this day
        intervals = []
        # TODO: Catch the RiseSetError for circumpolar stars
        # TODO: Return either a complete or empty interval, as appropriate
        # TODO: This requires either introspecting state in the error, or
        # TODO: calling an is_circumpolar method before this calculation

        if target == 'sun':
            transits, rises, sets = calc_sunrise_set(self.site, dt, self.twilight)
        elif target == 'moon':
            transits, rises, sets = calc_planet_rise_set(self.site, dt, MOON_REFRACTION, 'moon')
        else:
            # Test for circumpolarity
            if star.is_always_up(dt):
                # Return a full interval over the entire day
                intervals.append((dt, dt + ONE_DAY))
                return intervals

            # Catch target never rising
            try:
                effective_horizon_angle = Angle(degrees=effective_horizon)
                transits, rises, sets = calc_rise_set(target, self.site,
                                                      dt, effective_horizon_angle)
            except RiseSetError:
                return intervals


        _log.debug("latitude: %s",     self.site['latitude'].in_degrees())
        _log.debug("longitude: %s",    self.site['longitude'].in_degrees())
        _log.debug("twilight: %s",     self.twilight)
        _log.debug("dt: %s",           dt)
        _log.debug("rise: %s (%s)",    rises, dt + rises)
        _log.debug("transit: %s (%s)", transits, dt + transits)
        _log.debug("set: %s (%s)",     sets, dt + sets)

#        import ipdb; ipdb.set_trace()

        # Case 1: Overlapping start of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #         |       x                                     |
        #         |    x     x                                  |
        #         | x           x                               | x
        #        x|                x                           x|
        #     x   |                   x                     x   |
        #   Rise  0hr  Transit       Set                  Rise 24hr
        if (rises > transits) and (sets > transits):

            # Store the first interval - start of day until target set
            absolute_set = sets + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rises + dt
            intervals.append((absolute_rise, dt + ONE_DAY))


        # Case 2: Rise, set and transit all fall within the day, in order
        # Target rises today, transits, and sets before the day ends
        #         |                      x                     |
        #         |                   x     x                  |
        #         |                x           x               |
        #         |             x                 x            |
        #         |          x                       x         |
        #         |       x                             x      |
        #         0hr   Rise           Transit          Set   24hr
        elif (rises < transits) and (sets > transits):
            # Only one interval - rise until target set
            absolute_rise = rises + dt
            absolute_set  = sets + dt
            intervals.append((absolute_rise, absolute_set))


        # Case 3: Overlapping end of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #                 x    |                                        x    |
        #              x     x |                                     x     x |
        #           x          |x                                 x          |x
        #        x             |   x                           x             |   x
        #     x                |      x                     x                |
        #   Rise       Transit 0hr   Set                  Rise      Transit 24hr
        elif (rises < transits) and (sets < transits):
            # Same code as case 1!
            # Store the first interval - start of day until target set
            absolute_set = sets + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rises + dt
            intervals.append((absolute_rise, dt + ONE_DAY))


        return intervals


    def __repr__(self):
        repr_dict = copy.deepcopy(self.__dict__)
        del(repr_dict['dark_intervals'])

        sorted_str_dict = "{" + ", ".join("%s: %s" % (key, self.__dict__[key])
                             for key in sorted(self.__dict__)) + "}"
        return "Visibility (%s)" % sorted_str_dict

    def __key(self):
        return (self.site, self.start_date, self.end_date, self.horizon, self.twilight)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())
