import glob

files = glob.glob('./Patches/*.cpp')
print(files)


#with open(str(files), 'r') as f:
for i in range (len(files)):
  print (str(files[i]))
  f=open(str(files[i]), 'r')
  lines = f.readlines()
  print (lines)
