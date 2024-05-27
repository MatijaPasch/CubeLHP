# in this program, we check which normalised configurations were witnessed and store the remaining configurations in missing_configs.txt

import json

# program entry point
# we load the normalised configurations from normalised_configs.txt
# we load the witnessed configurations from witnessed_configs.txt
# we store the normalised configurations that were not witnessed in missing_configs.txt

with open('normalised_configs.txt','r') as file:
  normalised_configs=json.load(file)

with open('witnessed_configs.txt','r') as file:
  witnessed_configs=json.load(file)

missing_configs=[]
for i in range(len(normalised_configs)):
  if normalised_configs[i] not in witnessed_configs:
    missing_configs.append(normalised_configs[i])
with open('missing_configs.txt','w') as file:
  json.dump(missing_configs,file)
