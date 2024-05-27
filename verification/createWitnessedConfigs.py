# in this program, we verify that all elements in witnesses.txt are loose paths of maximal length minus 1
# for the elements which are such loose paths (meaning all of them), we write the corresponding configurations, meaning [start, end, skipped 1, skipped 2], to the file "witnessed_configs.txt"
#
# We use the following encoding of loose paths
# each loose path is given by 81 points, namely all points in the cube hypergraph
# the first hyperedge in the path is given by points 1, 2 and 3.
# the second hyperedge in the path is given by points 3, 4 and 5.
# The i-th hyperedge in the path is given by points 2(i-1)+1, 2(i-1)+2 and 2(i-1)+3 for i in {1,...,39}
# Point 80 is the first skipped point, point 81 is the second skipped point.

import json

# we test if the list lp is a loose path
# first, we check that it covers 81 elements (which are the points)
# then, we check that each point of the cube hypergraph is contained in lp
#       this means that lp contains 81 distinct points which are exactly the points of the cube hypergraph
# then, for each hyperedge we check that on exactly three coordinates the values of all three points coincide
#       and on exactly one coordinate the three points cover all three values
#       Implementation:
#         Let x=(x1,x2,x3,x4), y=(y1,y2,y3,y4), z=(z1,z2,z3,z4) be the three points.
#         Let ni=|{xi,yi,zi}| for i in {1,...,4} be the number of values attained
#         Let f=(|{i:ni=n}|)_n over n in {0,1,2,3} be the corresponding frequencies
#         Then we need f=(0,3,0,1)
# If these conditions are met, then lp indeed encodes a loose path
def is_loose_path(lp):
  if len(lp) != 81:
    return False
  cube_points=[[x0,x1,x2,x3] for x0 in range(3) for x1 in range(3) for x2 in range(3) for x3 in range(3)]
  present=[False for _ in range(81)]
  for pt in lp:
    if pt in cube_points:
      present[cube_points.index(pt)]=True
    else:
      return False
  if not all(present):
    return False
  ind =0
  for _ in range(39):
    pt1=lp[ind]
    pt2=lp[ind+1]
    pt3=lp[ind+2]
    val_cnt_freq=[0,0,0,0]
    for i in range(4):
      val_cnt=len(set([pt1[i],pt2[i],pt3[i]]))
      val_cnt_freq[val_cnt]=val_cnt_freq[val_cnt]+1
    if val_cnt_freq != [0,3,0,1]:
      return False
    ind=ind+2
  return True

# program entry point
# we load the witnesses, then for each witness we store the configuration (start, end, skipped 1, skipped 2) if the witness encodes a loose path (which it always does)
# then we store the witnessed configurations in witnessed_configs.txt

with open('witnesses.txt','r') as file:
  loose_paths=json.load(file)

cnt=len(loose_paths)
witnessed_configs=[]
for i in range(cnt):
  if is_loose_path(loose_paths[i]):
    witnessed_configs.append([loose_paths[i][0],loose_paths[i][78],loose_paths[i][79],loose_paths[i][80]])

with open('witnessed_configs.txt','w') as file:
  json.dump(witnessed_configs,file)
