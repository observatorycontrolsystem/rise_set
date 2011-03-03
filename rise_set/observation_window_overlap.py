# slampoud@cs.ucsb.edu 3/1/2011

from operator import attrgetter, methodcaller

class ObservationEndpoint:
  def __init__(self, time, amistart, seqid):
	self.time=time
        self.amistart=amistart
        self.seqid=seqid
  def __repr__(self):
      return repr((self.time, self.amistart, self.seqid))
  def getTime(self):
        return self.time
  def amIStart(self):
        return self.amistart
  def getSeqID(self):
        return self.seqid

# populate two lists of ObservationEndpoint objects with starts and ends

seq1=[]
seq2=[]
seq1.append(ObservationEndpoint(1, 1, 0)) # start
seq1.append(ObservationEndpoint(2, 0, 0)) # end

seq2.append(ObservationEndpoint(1, 1, 1)) # start
seq2.append(ObservationEndpoint(1.5, 0, 1)) # end
seq2.append(ObservationEndpoint(1.6, 1, 1)) # start
seq2.append(ObservationEndpoint(3, 0, 1)) # end

# concatenate sequences and sort by time 
# tie breaking: end and start w/ same time => end before start
#               2 ends or 2 starts w/ same time => order doesn't matter
# note: I don't know how to do double criterion sorting with methodcaller, which 
# would be the preferable way, so I'm using attrgetter
bothseq=sorted(seq1+seq2, key=attrgetter('time', 'amistart'))

overlaps=[]
flags=[0,0]
for oe in bothseq:
   if oe.amIStart():
      flags[oe.getSeqID()]+=1 # raise the flag
      if flags[(oe.getSeqID()+1)%2]: # detect start of an overlap
         overlaps.append(ObservationEndpoint(oe.getTime(), 1, 0)) # add start
   else:
      flags[oe.getSeqID()]-=1 # lower the flag
      if flags[(oe.getSeqID()+1)%2]: # detect end of an overlap
         overlaps.append(ObservationEndpoint(oe.getTime(), 0, 0)) # add end

print overlaps
