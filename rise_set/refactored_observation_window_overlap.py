# slampoud@cs.ucsb.edu 3/1/2011

from operator import attrgetter, methodcaller

class ObservationEndpoint:
  def __init__(self, time, is_start, seqid):
        self.time     = time
        self.is_start = is_start
        self.seqid    = seqid

  def __repr__(self):
      return repr( (self.time, self.is_start, self.seqid) )


# populate two lists of ObservationEndpoint objects with starts and ends

seq1 = []
seq2 = []
seq1.append(ObservationEndpoint(1, True, 0)) # start
seq1.append(ObservationEndpoint(2, False, 0)) # end

seq2.append(ObservationEndpoint(1, True, 1)) # start
seq2.append(ObservationEndpoint(1.5, False, 1)) # end
seq2.append(ObservationEndpoint(1.6, True, 1)) # start
seq2.append(ObservationEndpoint(3, False, 1)) # end

# concatenate sequences and sort by time
# tie breaking: end and start w/ same time => end before start
#               2 ends or 2 starts w/ same time => order doesn't matter
# note: I don't know how to do double criterion sorting with methodcaller, which 
# would be the preferable way, so I'm using attrgetter


def find_overlaps(seq1, seq2):
    bothseq = sorted(seq1+seq2, key=attrgetter('time', 'is_start'))

    overlaps = []
    flags    = [0,0]
    for oe in bothseq:
       if oe.is_start:
          flags[oe.seqid] += 1 # raise the flag
          if flags[(oe.seqid+1)%2]: # detect start of an overlap
             overlaps.append(ObservationEndpoint(oe.time, True, 0)) # add start
       else:
          flags[oe.seqid]-=1 # lower the flag
          if flags[(oe.seqid+1)%2]: # detect end of an overlap
             overlaps.append(ObservationEndpoint(oe.time, False, 0)) # add end

    return overlaps


print "Method 1"
overlaps = find_overlaps(seq1, seq2)
print overlaps


def windows_intersect(start1, end1, start2, end2):
    if end2 <= start1 or start2 >= end1:
        return False

    return  start1 <= end2 and end1 >= start2



def find_overlaps2(seq1, seq2):
    overlaps = ()

    return overlaps


seq3 = (
         (1, 2),
       )

seq4 = (
         (1, 1.5),
         (1.6, 3),
       )


print "\nAlternative method"
overlaps = find_overlaps2(seq3, seq4)
print overlaps
