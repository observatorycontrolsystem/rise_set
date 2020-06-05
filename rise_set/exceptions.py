""" Exceptions used throughout rise-set package """
class InvalidAngleError(Exception):

    """ 
    Exception that is raised when out-of-range angles are provided to the
    Angle class.
    """
    pass


class AngleConfigError(Exception):

    """
    Exception that is raised when invalid constructor arguments are given
    to the Angle class.
    """
    pass

class InvalidDateTimeError(Exception):

    """ Exception that is raised when an invalid date is encountered."""
    pass


class RiseSetError(Exception):

    """ Exception that is raised when a target either never rises or never sets."""
    pass


class IncompleteTargetError(Exception):

    """ Exception that is raised when a target is missing a key value (RA, Dec, etc.)."""
    pass


class RatesConfigError(Exception):

    """ 
    Exception that is raised when invalid constructor arguments
    are provided to the ProperMotion class.
    """
    pass

class MovingViolation(Exception):

    """ Exception that is raised when moving object errors are encountered."""
    pass

class InvalidHourAngleLimit(Exception):
    
    """ Exception that is raised when an invalid hour angle limit is given."""
    pass
