
Matrix2=np.array(Matrix)#Numpy matrix is very flexible. ie a=Weights[295,0:400]  
Weights=np.array(Weights)

#Do have to declare some double variables in python. 

#Perhaps its only lists that don't need declaring. Their declaration is implied by their use syntax however for doubles. It is not possible that the context of use implies declaration. Its too ambigous.
column_index=0
row_index=0
old_row=0
old_column=0

for i in xrange(0,int(Matrix2.shape[1])):
 for j in xrange(0,int(Weights.shape[0])):
  if sum(Weights[j,:]) > old_row: # old :
   old_row = sum(Matrix2[j,:])
   row_index=j
  if sum(Matrix2[:,i]) > old_column: # old :
   old_column = sum(Matrix2[:,i])
   column_index=i
     
     
outdegree=column_index
indegree=row_index

column_index=0
row_index=0
old_row=0
old_column=0
weightdegree=0

for i in xrange(0,int(Weights.shape[1])):
 for j in xrange(0,int(Weights.shape[0])):
  if sum(Weights[j,:]) > old_row: # old :
   old_row = sum(Weights[j,:])
   row_index=j
  if sum(Weights[:,i]) > old_column: # old :
   old_column = sum(Weights[:,i])
   column_index=i
     
weightdegree=column_index

#weights are sent via netcon. at the source. Does not make sense its a matrix.
#Its only a column vector.

#print column_index, "high in degree cell"
#print row_index, "high out degree cell"
