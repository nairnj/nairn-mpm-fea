import glob

files = glob.glob('./Global_Quantities/*.cpp')
# print(files)


#with open(str(files), 'r') as f:
for i in range (len(files)):
  # print (str(files[i]))
  f=open(str(files[i]), 'r+')
  fs=open(str(files[i])+"_", 'w')
  lines = f.readlines()
  for l in range (len(lines)):
    if (lines[l].find("stdafx.h") != -1):
      print ("found \n" + str(files[i]))
    #fs.write
      fs.write("#if defined ( _MSC_VER) || defined (__APPLE__)\n")
      fs.write (lines[l])
      fs.write ("#endif\n")
      
      print ("#if defined ( _MSC_VER) || defined (__APPLE__)\n")
      print (lines[l])
      print ("#endif")
    else:
      fs.write (lines[l])