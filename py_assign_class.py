
#for x in range(0, 3):
for i in range(0,len(name_list)-1):
  #print str(name_list[i])
  if "pyramid" not in str(name_list[i]) :continue 
  #h.forsec{ h.cells(i).all{ uninsert hh }
  #h.cells(i).ca1pyrHIPP()
  
  print "class pyr", i
  
for i in range(0,len(name_list)-1):
  print name_list[i]  
  
  
for i in range(0,len(name_list)-1):  
  if "spindle" not in str(name_list[i]) :continue 
  print "class spindle", i
  
##Fast Spiking Inhibitory.  
for i in range(0,len(name_list)-1):  
  if "basket" not in str(name_list[i]) :continue 
  print "basket", i  
