#!/usr/bin/python

import os
import sys

if len(sys.argv) < 3:
    print "Usage: %s <moviefile_basename> <fps>"%(sys.argv[0])
    sys.exit(1)

fps=int(sys.argv[2])
movie_filename=sys.argv[1]
try:
    image_prefix=sys.argv[3]
except:
    image_prefix=""

cmd="mencoder mf://%s*.png -mf fps=%d -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=1600:vhq -ffourcc XVID -o %s"%(image_prefix,fps,movie_filename+".avi")
os.system(cmd)
