# in this program, we determine all normalised configurations and write them to the file normalised_configs.txt

import json

# https://math.stackexchange.com/q/4876838/1096797
# the row "row" of the configuration "config" is relabeled with respect to the order of first occurrence
# enum is the enumeration x of X={config[i][row]:i} with respect to the order of first occurrence
def relabel_by_first_occurrence(config,row):
  enum=[config[0][row]]
  for i in range(1,4):
    if config[i][row] not in enum:
      enum.append(config[i][row])
  for i in range(4):
    config[i][row]=enum.index(config[i][row])

# relabels each single row of the configuration config
def relabel_values_per_coordinate(config):
  for i in range(4):
    relabel_by_first_occurrence(config,i)

# sort the rows of the configuration "config" lexicographically
def sort_rows(config):
  sorted_transposed=sorted([[config[j][i] for j in range(4)] for i in range(4)],key=lambda x: (x[0],x[1],x[2],x[3]))
  return [[sorted_transposed[j][i] for j in range(4)] for i in range(4)]

# relabel the rows of the configuration "config" independently and then sort the rows lexicographically
# these are the first two steps of the normalisation of "config"
# the last step where these first two steps are applied to "config" with the skipped points switched is excluded here
def normalise_excl_skip_switch(config):
  relabel_values_per_coordinate(config)
  return sort_rows(config)

# consider config1 and config2 in matrix notation where the points are the columns
# then we consider the order on the 16-dim points given by the juxtaposition of the rows
# we return -1 if config1 is smaller, 0 if they are equal and 1 otherwise
def compare_configurations(config1,config2):
  for i in range(4):
    for j in range(4):
      if config1[j][i]<config2[j][i]:
        return -1
      if config1[j][i]>config2[j][i]:
        return 1
  return 0

# the configuration "config" is normalised
def normalise(config):
  config_skip_switched=[config[0],config[1],config[3],config[2]]
  config_normalised=normalise_excl_skip_switch(config)
  config_skip_switched_normalised=normalise_excl_skip_switch(config_skip_switched)
  if compare_configurations(config_normalised,config_skip_switched_normalised)>0:
    return config_skip_switched_normalised
  return config_normalised

# checks if the four points in the configuration "config" are pairwise distinct
def points_distinct(config):
  for i in range(1,4):
    for j in range(i):
      if config[i]==config[j]:
        return False
  return True

# generate all normalised configurations
# we achieve this by cycling through all configurations (with distinct points)
# then we normalise them, yielding all normalised configurations
def generate_normalised_configs():
  normalised_configs=[]
  all_points=[[x1,x2,x3,x4] for x1 in range(3) for x2 in range(3) for x3 in range(3) for x4 in range(3)]
  point_cnt= 3 ** 4
  config_cnt=int( (point_cnt ** 3) *(point_cnt-1) / 2 )
  config_ind=0
  for start_point in all_points:
    for end_point in all_points:
      for ind_s1 in range(1,point_cnt):
        for ind_s2 in range(ind_s1):
          if config_ind % 100000 == 0:
            print(f"\r{round(config_ind/config_cnt*100,2)} % complete",end=' ',flush=True)
          config_ind=config_ind+1
          skipped_point1=all_points[ind_s1]
          skipped_point2=all_points[ind_s2]
          config=[list(start_point),list(end_point),list(skipped_point1),list(skipped_point2)]
          if points_distinct(config):
            normalised_config=normalise(config)
            if normalised_config not in normalised_configs:
              normalised_configs.append(normalised_config)
  return normalised_configs


# program entry point
# we generate all normalised configurations and write them to the file normalised_configs.txt
normalised_configs=generate_normalised_configs()
print("\r100.00 % complete")
print(f"{len(normalised_configs)} normalised configurations found")
with open('normalised_configs.txt','w') as file:
  json.dump(normalised_configs,file)