import time
from stinger import *
from random import *

nv = 262144
ne = nv * 8

sv = int64Array(ne)
ev = int64Array(ne)
w = int64Array(ne)

print "generating edges"

for i in range(ne):
  sv[i] = randint(0, nv-1)
  ev[i] = randint(0, nv-1)
  w[i] = 1

print "creating stinger"

s = edge_list_to_stinger(nv, ne, sv, ev, w, None, None, 0)

print "new iterator"
it = stinger_iterator_new(s)
types = int64Array(1)
types[0] = 0

print "setup iterator"
stinger_iterator_edge_type_filter(types, 1, it)

print "calculating components"
changed = 0
components = range(nv)
start = time.clock()
while(1):
  changed = 0
  while(stinger_iterator_next(it)):
    if(components[it.source] < components[it.dest]):
      components[it.dest] = components[it.source]
      changed = 1
  if(not changed):
    break
  for i in range(nv):
    while(components[i] != components[components[i]]):
      components[i] = components[components[i]]
stop = time.clock()

print "done. time = " + str(stop - start)

print "counting components"

count = 0
for i in range(nv):
  if(components[i] == i):
    count += 1

print count
