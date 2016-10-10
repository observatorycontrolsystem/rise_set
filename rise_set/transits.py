#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import next
from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
import ast
from datetime import datetime,timedelta
from visibility import Visibility
from sky_coordinates import RightAscension, Declination
from angle import Angle
import sys, argparse
from math import sqrt
from reqdb.utils.duration import calculate_duration
from reqdb.requests import Request, UserRequest
from reqdb.client   import SchedulerClient

# breaks up long intervals into discrete chunks that are schedulable
# by site with a minimum overlap of MINOVERLAP minutes and a maximum
# overlap of MAXOVERLAP minutes. Negative values of MINOVERLAP imply
# maximum gaps are acceptable. By default, it will schedule everything possible.

def distribute(interval, target, minoverlap=-1440, maxoverlap=1440, onepersite=True):

    sitenames = []
    intervals = []
    titles = []
    telescopes = file_to_dicts('/home/jeastman/lcogt/scheduler/trunk/test/telescopes.dat')
    for telescope in telescopes:
        if telescope['status'] == 'online':
            site =  {
                'name'     : telescope['name'].split(".")[2],
                'latitude' : Angle(degrees=telescope['latitude']),
                'longitude': Angle(degrees=telescope['longitude']),
                'horizon'  : telescope['horizon'],
                'ha_limit_neg' : telescope['ha_limit_neg'],
                'ha_limit_pos' : telescope['ha_limit_pos'],
                }            

            if not site['name'] in sitenames or not onepersite:
                sitenames.append(site['name'])
                print(site['name'])
                for obs_start,obs_end in interval:
                    observability = Visibility(
                        site=site,
                        start_date=obs_start,
                        end_date=obs_end,
                        horizon=site['horizon'],
                        twilight='nautical',
                        ha_limit_neg=site['ha_limit_neg'],
                        ha_limit_pos=site['ha_limit_pos'],
                        )
                    site_intervals = observability.get_observable_intervals(target, moon_distance=Angle(degrees=0))
                    for start,stop in site_intervals:
                        intervals.append((start,stop))
                        titles.append(target['name'] + ' (' + site['name'] + ')')

    return intervals, titles

# Read a file of code into a list (for telescopes.dat)
def file_to_dicts(filename):
    fh = open(filename, 'r')
    data = fh.read()
    return ast.literal_eval(data)

# read in the exoplanets.org csv file
def readexo():
    import csv
    f = open('exoplanets.csv','rb')
    reader = csv.reader(f)
    headers = next(reader)
    planets = {}
    for h in headers:
        planets[h] = []
    for row in reader:
        for h,v in zip(headers,row):
            planets[h].append(v)
    return planets


# get the number of exposures we can take during a given window for a single molecule
def get_nexp(window, molecule):

    diff = window['end'] - window['start']
    window_duration = old_div(diff.microseconds,1e6) + diff.seconds + diff.days*86400.0
    molecule['exposure_count'] = int(old_div(window_duration,molecule['exposure_time']))

    while True:
        molecule_duration = calculate_duration([molecule])
        if molecule_duration.seconds <= window_duration:
            break
        molecule['exposure_count'] -= 1


# Find the exoplanet parameters required to calculate the ephemeris from exoplanets.org
def get_ephem(name):

    # The name is specified, use values from exoplanets.org
    planets = readexo()
                    
    for i in range(len(planets['NAME'])):
        if planets['NAME'][i] == opt.name:
            # we should ensure these exist prior to using them
            per = float(planets['PER'][i])
            uper = float(planets['UPER'][i])
            tt = float(planets['TT'][i])
            utt = float(planets['UTT'][i])
            t14 = float(planets['T14'][i])
            ut14 = float(planets['UT14'][i])
            ra = RightAscension(degrees=float(planets['RA'][i])*15.0)
            dec = Declination(degrees=float(planets['DEC'][i]))
            vmag = float(planets['V'][i])

            target =  {
                'name'              : opt.name,
                'ra'                : ra,
                'dec'               : dec,
                'epoch'             : 2000,
                }
           
            return per, uper, tt, utt, t14, ut14, ra, dec, vmag
        
    # not found, return an error
    return -1

# calculates the transit intervals for all sites
def get_all_transit_intervals(target, ephemeris, start_date, end_date, onepersite=False, verbose=False, vverbose=False,exptime=0):
     
    intervals = []
    titles=[]
    sitenames = []
    telescopes = file_to_dicts('/home/jeastman/lcogt/scheduler/transits/rise_set-0.2.10/telescopes.dat')
    for telescope in telescopes:
        if telescope['status'] == 'online':
            site =  {
                'name'     : telescope['name'].split(".")[2],
                'latitude' : Angle(degrees=telescope['latitude']),
                'longitude': Angle(degrees=telescope['longitude']),
                'horizon'  : telescope['horizon'],
                'ha_limit_neg' : telescope['ha_limit_neg'],
                'ha_limit_pos' : telescope['ha_limit_pos'],
                }
            if not site['name'] in sitenames or not onepersite:
                sitenames.append(site['name'])

                site_intervals,site_titles = get_transit_intervals(target, ephemeris, site, start_date, end_date, verbose=verbose, vverbose=vverbose,exptime=exptime)

#                print site_intervals, site_titles, start_date, end_date, target
                for start,stop in site_intervals:
                    intervals.append((start,stop))
                for title in site_titles:
                    titles.append(title)

    return intervals,titles
    

def get_transit_intervals(target, ephemeris, site, start_date, end_date, verbose=False, vverbose=False,exptime=0):

    # period in days (float)
    # tt is a datetime object of the center of the transit (BJD_TDB)
    # t14 is the duration of transit, in days
    # target is a target object that must include ra, dec
    # start_date is a datetime object of the start of the window (UTC)
    # end_date is a datetime object of the end of the window (UTC)
    
    # calculate the epoch numbers that correspond to the observing window    
    startdiff = start_date-ephemeris['tt']
    startepoch = int(old_div((old_div(startdiff.microseconds,864e10) + old_div(startdiff.seconds,86400) + startdiff.days),ephemeris['period']))
    enddiff = end_date-ephemeris['tt']
    stopepoch = int(old_div((old_div(enddiff.microseconds,864e10) + old_div(enddiff.seconds,86400) + enddiff.days),ephemeris['period']))+1

    sitestr = site['name']
    visible_transits = []
    titles = []
    for epoch in range(startepoch, stopepoch+1):
        
        ## Uncertainty in the Transit Time of the Nth transit
        uttn = sqrt(ephemeris['utt']**2 + (epoch*ephemeris['uperiod'])**2 + ephemeris['ut14']**2)
        
        ## this defines the ideal window, including 3-sigma uncertainty
        contact1 = ephemeris['tt'] + timedelta(days=ephemeris['period']*epoch - 0.5*ephemeris['t14']) #- 3*uttn
        contact4 = ephemeris['tt'] + timedelta(days=ephemeris['period']*epoch + 0.5*ephemeris['t14']) #+ 3*uttn
        event_start = contact1 - timedelta(days=(ephemeris['maxpre' ])*ephemeris['t14'] + 3.0*uttn)
        event_stop  = contact4 + timedelta(days=(ephemeris['maxpost'])*ephemeris['t14'] + 3.0*uttn)
#        print (ephemeris['maxpre' ]),ephemeris['t14'], 3.0*uttn, event_stop-event_start, event_start,event_stop

        ## convert event_start, event_stop, contact1, contact4 from BJD_TDB to JD_UTC
        # ignore for now; this will introduce a systematic error of up to ~9 minutes (depending on target and time of year)
        
        # determine when the object is observable from this telescope
        observability = Visibility(
            site=site,
            start_date=max(event_start,start_date),
            end_date=min(event_stop,end_date),
            horizon=site['horizon'],
            twilight='nautical',
            ha_limit_neg=site['ha_limit_neg'],
            ha_limit_pos=site['ha_limit_pos'],
            )
        observable_intervals = observability.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        for visible in observable_intervals:
    
            # the maximum we could observe it
            obs_start = max(event_start,visible[0])
            obs_end = min(event_stop,visible[1])
            
            # how much of the transit is visible?
            transitvis = min(contact4,visible[1])-max(contact1,visible[0])
            transitfrac = max(round(old_div((old_div(transitvis.microseconds,864e10) + old_div(transitvis.seconds,86400.0) + transitvis.days),ephemeris['t14']),3),0)

            # how much prior to ingress is visible?
            previs = min(contact1,visible[1]) - max(event_start,visible[0])
            prefrac = round(max(old_div((old_div(previs.microseconds,864e10) + old_div(previs.seconds,86400.0) + previs.days),ephemeris['t14']),0),3)

            # how much after egress is visible?
            postvis = min(event_stop,visible[1])-max(contact4,visible[0])
            postfrac = round(max(old_div((old_div(postvis.microseconds,864e10) + old_div(postvis.seconds,86400.0) + postvis.days),ephemeris['t14']),0),3)

            # how much out of transit total is visible?
            ootfrac = prefrac + postfrac

#            print verbose

            # does it meet the criteria?
            if transitfrac >= ephemeris['minevent'] and prefrac >= ephemeris['minpre'] and \
                    postfrac >= ephemeris['minpost'] and ootfrac >= ephemeris['minbaseline']:
                visible_transits.append((obs_start,obs_end))
                titles.append(target['name'] + ' ' + site['name'] + ' ' + str(epoch) + ' ' + str(transitfrac) + ' ' + str(ootfrac))
                if verbose or vverbose:
                    window = {
                        'start' : obs_start + timedelta(minutes=5),
                        'end' : obs_end - timedelta(minutes=5)
                        }    
                    molecule = {
                        # Required fields
                        'exposure_time'   : exptime,        # Exposure time, in secs
                        'exposure_count'  : 0,              # The number of consecutive exposures
                        'filter'          : '',           # The generic filter name
                        # Optional fields. Defaults are as below.
                        'type'            : 'EXPOSE',       # The type of the molecule
                        'ag_name'         : 'SCICAM',       # '' to let it resolve; 'SCICAM' for self-guiding
                        'ag_mode'         : 'Optional',
                        'instrument_name' : 'SCICAM',       # This resolves to the main science camera on the scheduled resource
                        'bin_x'           : 2,              # Your binning choice. Right now these need to be the same.
                        'bin_y'           : 2,
                        }
#                    get_nexp(window, molecule)
                    print(site['name'], epoch, obs_start, obs_end, contact1, contact4, transitfrac, ootfrac, uttn*1440.0)
            else:
                if vverbose:
                    print(site['name'], epoch, obs_start, obs_end, contact1, contact4, transitfrac, ootfrac, molecule['exposure_count'], uttn*1440.0, "*")

                
    return visible_transits,titles

def get_exptime(vmag):
    defocus = 0.0
    m0 = 12.0
    exptime = 25*10.0**(old_div((vmag-m0),2.5))
    exptime = min(round(exptime,0),300)
    while exptime < 30:
        exptime += 10.0
        defocus += 1.0
        
    return exptime, defocus

def schedule_intervals(intervals, target, band, vmag, propid, email, titles, verbose=False, exptime=None, defocus=None):
    constraints = {'max_airmass' : 2.0}

    location = {'telescope_class' : '1m0'}

    proposal = {
        'proposal_id'   : propid,
        'user_id'       : email,
        }


    if exptime == None:
        exptime,defocus = get_exptime(vmag)
    elif defocus==None:
        defocus = 0

    # if defocused, guide with the guider
    # otherwise, use the science camera
    if defocus > 0:
        guider = ''
    else:
        guider = 'SCICAM'

    molecule = {
        # Required fields
        'exposure_time'   : exptime,        # Exposure time, in secs
        'exposure_count'  : 0,              # The number of consecutive exposures
        'filter'          : band,           # The generic filter name
        # Optional fields. Defaults are as below.
        'type'            : 'EXPOSE',       # The type of the molecule
        'ag_name'         : 'SCICAM',       # '' to let it resolve; 'SCICAM' for self-guiding
        'ag_mode'         : 'Optional',
        'instrument_name' : guider,         # This resolves to the main science camera on the scheduled resource
        'bin_x'           : 2,              # Your binning choice. Right now these need to be the same.
        'bin_y'           : 2,
        'defocus'         : defocus         # Mechanism movement of M2, or how much focal plane has moved (mm)
        }
    
    # kludge because target changes definitions between Visibility and set_target
    target['ra'] = target['ra'].degrees
    target['dec'] = target['dec'].degrees
 
    i=0
    for obs_start,obs_end in intervals:

        # make the window a little shorter to calculate the number of exposures
        window = {
            'start' : obs_start + timedelta(minutes=2),
            'end' : obs_end - timedelta(minutes=2)
            }
        get_nexp(window, molecule)
        molecule['exposure_count'] -= 2

        req = Request()
        req.set_location(location)
        req.set_target(target)

        # this is the actual window
        window = {
            'start' : obs_start,
            'end' : obs_end
            }    
        req.add_window(window)

        req.add_molecule(molecule)
        req.set_constraints(constraints)
        
        ur = UserRequest(group_id=titles[i])
        ur.add_request(req)
        ur.operator = 'single'
        ur.set_proposal(proposal)

        if verbose:
            print(molecule['exposure_time'], molecule['exposure_count'], molecule['defocus'], vmag, band, titles[i])


        # You're done! Send the complete User Request to the DB for scheduling
        client        = SchedulerClient('http://scheduler-dev.lco.gtn/requestdb/')
        response_data = client.submit(ur, keep_record=True, debug=True)
        client.print_submit_response()
        i = i+1

    return

if __name__ == '__main__':

    # EXAMPLES:
    # see visible transits of WASP-12b in the next month
    # python transits2.py --name=WASP-12b

    # see visible transits of WASP-12b between 1/1/2013 and 1/5/2013
    # python transits2.py --name=WASP-12b --startwin=20130101 --stopwin=20130105

    # Schedule transits of WASP-12b between 1/1/2013 and 1/5/2013
    # python transits2.py --name=WASP-12b --startwin=20130101 --stopwin=20130105 --schedule --band=ip --email=email@lcogt.net --propid=XXXXX-###

    # Schedule transits of NEW_OBJECT between 1/1/2013 and 1/5/2013 using your own ephemeris and target info
    # python transits2.py --name=NEW_OBJECT --startwin=20130101 --stopwin=20130105 --period=1.09142245 --t14=0.122 --tt=2454508.976820 --ra=97.6366416613 --dec=29.672277778 --schedule --band=ip --email=email@lcogt.net --propid=XXXXX-###
    

    # define the arguments and defaults
    today = datetime.strftime(datetime.today()+timedelta(days=1), format='%Y%m%d')
    month = datetime.strftime(datetime.today()+timedelta(days=30), format='%Y%m%d')
    parser = argparse.ArgumentParser(description='Schedule transit observations on the LCOGT network')
    parser.add_argument('--name'         ,dest='name'         ,action='store'     ,type=str                ,help='Planet name from exoplanets.org from which to get ra, dec, period, uperiod, tt, utt, t14, ut14. Specifying those separately will override the exoplanets.org values. Spaces are ignored -- e.g., "WASP-12b" is equivalent to "WASP-12 b"')
    parser.add_argument('--period'       ,dest='per'          ,action='store'     ,type=float              ,help='Period of the planet (days).')
    parser.add_argument('--uperiod'      ,dest='uper'         ,action='store'     ,type=float              ,help='Uncertainty in the period (days)')
    parser.add_argument('--tt'           ,dest='tt'           ,action='store'     ,type=float              ,help='Transit time of the 0th epoch (BJD_TDB)')
    parser.add_argument('--utt'          ,dest='utt'          ,action='store'     ,type=float              ,help='Uncertainty in the transit time (days)')
    parser.add_argument('--t14'          ,dest='t14'          ,action='store'     ,type=float              ,help='Duration of the transit, from 1st to 4th contact (days)')
    parser.add_argument('--ut14'         ,dest='ut14'         ,action='store'     ,type=float              ,help='Uncertainty in the duration (days)')
    parser.add_argument('--ra'           ,dest='ra'           ,action='store'     ,type=float              ,help='Right ascension of the target (J2000 degrees)')
    parser.add_argument('--dec'          ,dest='dec'          ,action='store'     ,type=float              ,help='Declination of the target (J2000 degrees)')
    parser.add_argument('--mintransit'   ,dest='mintransit'   ,action='store'     ,type=float,default=1.0  ,help='Minimum fraction of the transit that must be visible to be scheduled')
    parser.add_argument('--minoot'       ,dest='minoot'       ,action='store'     ,type=float,default=0.5  ,help='Minimum fraction of the out of transit time that must be visible to be scheduled, in units of T14')
    parser.add_argument('--minpre'       ,dest='minpre'       ,action='store'     ,type=float,default=0.0  ,help='Minimum fraction of the out of transit before ingress that must be visible to be scheduled, in units of T14')
    parser.add_argument('--minpost'      ,dest='minpost'      ,action='store'     ,type=float,default=0.0  ,help='Minimum fraction of the out of transit after egress that must be visible to be scheduled, in units of T14')
    parser.add_argument('--maxpre'       ,dest='maxpre'       ,action='store'     ,type=float,default=0.5  ,help='Maximum fraction of the out of transit before ingress to attempt to schedule, in units of T14')
    parser.add_argument('--maxpost'      ,dest='maxpost'      ,action='store'     ,type=float,default=0.5  ,help='Maximum fraction of the out of transit after egress to attempt to schedule, in units of T14')
    parser.add_argument('--startwin'     ,dest='startwin'     ,action='store'     ,type=str  ,default=today,help='UT date of the beginning of the observing window (YYYYMMDD). Default is today.')
    parser.add_argument('--stopwin'      ,dest='stopwin'      ,action='store'     ,type=str  ,default=month,help='UT date of the end of the observing window (YYYYMMDD). Default is next month.')
    parser.add_argument('--schedule'     ,dest='schedule'     ,action='store_true'           ,default=False,help='Schedule these transits (requires BAND, PROPID, EMAIL, VMAG)')
    parser.add_argument('--band'         ,dest='band'         ,action='store'     ,type=str                ,help='The fitler to observe the transit')
    parser.add_argument('--propid'       ,dest='propid'       ,action='store'     ,type=str                ,help='The proposal id to associate with the scheduler')
    parser.add_argument('--email'        ,dest='email'        ,action='store'     ,type=str                ,help='The email associated with the proposal')
    parser.add_argument('--vmag'         ,dest='vmag'         ,action='store'     ,type=float              ,help='The absolute V band magnitude of the target')
    parser.add_argument('--verbose'      ,dest='verbose'      ,action='store_true'           ,default=False,help='Print extra information')
    parser.add_argument('--vverbose'     ,dest='vverbose'     ,action='store_true'           ,default=False,help='Print even more information (for debugging)')
    parser.add_argument('--onepersite'   ,dest='onepersite'   ,action='store_true'           ,default=False,help='Only schedule one transit at each site')
# these haven't been implemented yet, but are on my to-do list...
    parser.add_argument('--exptime'      ,dest='exptime'      ,action='store'     ,type=float              ,help='The desired exposure time (either this or VMAG required if SCHEDULE is set)')
    parser.add_argument('--defocus'      ,dest='defocus'      ,action='store'     ,type=float              ,help='The defocus, in mm. Only used if SCHEDULE is set. If not set, it is calculated based on VMAG')
#    parser.add_argument('--secondary'    ,dest='secondary'    ,action='store_true',           default=False,help='Calculate secondary eclipses, not primary transits')
#    parser.add_argument('--nschedule'    ,dest='nschedule'    ,action='store'     ,type=int                ,help='The maximum number of transits to schedule. Default is all.')
#    parser.add_argument('--schedulebest' ,dest='schedulebest' ,action='store_true'           ,default=False,help='Only schedule NSCHEDULE best (sorted by faction of transit visible, then by fraction of OOT visible). Default is soonest.')
#    parser.add_argument('--multiplesites',dest='multiplesites',action='store_true',           default=False,help='Allow observations at multiple sites to satisfy constraints')
#    parser.add_argument('--minoverlap'   ,dest='minoverlap'   ,action='store'     ,type=float,default=0.0  ,help='Minimum overlap between sites for each transit (hrs). Default is 0.')
#    parser.add_argument('--maxoverlap'   ,dest='maxoverlap'   ,action='store'     ,type=float,default=24.0 ,help='Maximum overlap between sites for each transit (hrs). Default is 24.')
#    parser.add_argument('--search'       ,dest='search'       ,action='store_true'           ,default=False,help='Schedule observations that cover the entire 3-sigma uncertainty window of a single transit (use to search for transits)')

    opt = parser.parse_args()

    per = None
    tt = None
    t14 = None
    ra = None
    dec = None
    
    # it's an error not to specify either the name or period, tt, t14, ra, and dec
    # we should check for that...

    # The name is specified, use values from exoplanets.org
    if opt.name != None:

        planets = readexo()
                    
        found = False
        for i in range(len(planets['NAME'])):
            if planets['NAME'][i].replace(' ','') == opt.name.replace(' ',''):

                # get the ephemeris values from exoplanets.org
                # if blank, set to None (to overwrite/check later)
                if planets['PER'][i] == '':
                    per = None
                else:
                    per = float(planets['PER'][i])
                if planets['UPER'][i] == '':
                    uper = 0.0
                else:
                    uper = float(planets['UPER'][i])
                if planets['TT'][i] == '':
                    tt = None
                else:                   
                    tt = float(planets['TT'][i])
                if planets['UTT'][i] == '':
                    utt = 0.0
                else:
                    utt = float(planets['UTT'][i])
                if planets['T14'][i] == '':
                    t14 = None
                else:
                    t14 = float(planets['T14'][i])
                if planets['UT14'][i] == '':
                    t14 = 0.0
                else:
                    ut14 = float(planets['UT14'][i])
                

                ra = RightAscension(degrees=float(planets['RA'][i])*15.0)
                dec = Declination(degrees=float(planets['DEC'][i]))
                vmag = float(planets['V'][i])
                print(vmag)
                found = True
                break
        if not found:
            print(opt.name + ' not found in exoplanets.org; using command line arguments')




    # Override exoplanets.org values with input values, if specified
    if opt.per != None:
        per = opt.per
    if opt.uper != None:
        uper = opt.uper
    elif uper == None:
        uper = 0.0

    if opt.vmag != None:
        vmag = opt.vmag
        
    if opt.tt != None:
        tt = opt.tt    
    if opt.utt != None:
        utt = opt.utt
    elif utt == None:
        utt = 0.0

    if opt.t14 != None:
        t14 = opt.t14
    if opt.ut14 != None:
        ut14 = opt.ut14
    else:
        ut14 = 0.0

    if opt.ra != None:
        ra = RightAscension(degrees=opt.ra)
    if opt.dec != None:
        dec = Declination(degrees=opt.dec)

#    if per==None or tt==None or t14==None or not ra<>None or not dec<>None:
#        print "ERROR: must specify either --name or per, tt, t14, ra, dec"
#        sys.exit(0)

    # Error checking
    if opt.schedule:
        if opt.band == None:
            print("Must specify --band if --schedule is set")
            sys.exit(0)
        if vmag == None:
            print("Must specify --vmag if --schedule is set")
            sys.exit(0)    
        if opt.propid == None:
            print("Must specify --propid if --schedule is set")
            sys.exit(0)
        if opt.email == None:
            print("Must specify --email if --schedule is set")
            sys.exit(0)

    print(vmag - 3)
            

    if per == None:
        print("Period not found; must specify --period")
        sys.exit(0)
    if tt == None:
        print("Transit time not found; must specify --tt")
        sys.exit(0)
    if t14 == None:
        print("Duration not found, must specify --t14")
        sys.exit(0)

#    if ra == None:
#        print "RA not found, must specify --ra"
#        sys.exit(0)
#    if len(dec) == 0:
#        print "Dec not found; must specify --dec"
#        sys.exit(0)

    if opt.vverbose:
        print(per, uper, tt, utt, t14, ut14, ra.degrees, dec.degrees)

    
    # define the target
    target = {
        'name'              : opt.name,
        'ra'                : ra,
        'dec'               : dec,
        'epoch'             : 2000,
        }

    # convert input dates to datetime objects
    start_date = datetime.strptime(opt.startwin,'%Y%m%d')
    end_date = datetime.strptime(opt.stopwin,'%Y%m%d')
    tt = datetime(1858,11,17,0,0,0) + timedelta(tt-2400000.5)

    ephemeris = {
        'tt'                : tt,
        'utt'               : utt,
        'period'            : per,
        'uperiod'           : uper,
        't14'               : t14,
        'ut14'              : ut14,
        'minevent'          : opt.mintransit,
        'minbaseline'       : opt.minoot,
        'minpre'            : opt.minpre,
        'minpost'           : opt.minpost,
        'maxpre'            : opt.maxpre,
        'maxpost'           : opt.maxpost,
        }

    print(ephemeris, target, ra.degrees)

    # get the transit intervals for each site
    intervals,titles = get_all_transit_intervals(target, ephemeris, start_date, end_date, onepersite=opt.onepersite, verbose=opt.verbose or not opt.schedule, vverbose=opt.vverbose, exptime=opt.exptime)
   
    # schedule the transits
    if opt.schedule:
        schedule_intervals(intervals, target, opt.band, vmag, opt.propid, opt.email, titles, verbose=opt.vverbose, exptime=opt.exptime, defocus=opt.defocus)
