 Matrix2=np.array(Matrix)#Numpy matrix is very flexible. ie a=Matrix2[295,0:400]  
 Weights=np.array(Weights)


 column_index=0
 row_index=0
 old_row=0
 old_column=0
 for i in xrange(0,int(Matrix2.shape[1])):
  for j in xrange(0,int(Matrix2.shape[0])):
   if sum(Matrix2[j,:]) > old_row: # old :
    old_row = sum(Matrix2[j,:])
    row_index=j
   if sum(Matrix2[:,i]) > old_column: # old :
    old_column = sum(Matrix2[:,i])
    column_index=i
     
     
outdegree=column_index
indegree=row_index     
#print column_index, "high in degree cell"
#print row_index, "high out degree cell"
