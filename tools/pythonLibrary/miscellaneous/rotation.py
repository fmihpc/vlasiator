import numpy as np

#I = np.array([[1,0,0],
#              [0,1,0],
#              [0,0,1]])
#B_vol = ..
#PTensorDiagonal = variables["PTensorDiagonal"]
#PTensorOffDiagonal = variables["PTensorOffDiagonal"]
#PTensor = [np.array([[PTensorDiagonal[i][0], PTensorOffDiagonal[i][2], PTensorOffDiagonal[i][1]],
#                     [PTensorOffDiagonal[i][2], PTensorDiagonal[i][1], PTensorOffDiagonal[i][0]],
#                     [PTensorOffDiagonal[i][1], PTensorOffDiagonal[i][0], PTensorDiagonal[i][2]]]) for i in xrange(len(PTensorDiagonal))]
#PTensor_rotated = []
#for i in xrange(len(B_vol)):
#   B_vol_iter = B_vol[i]
#   PTensor_iter = PTensor[i]
#   theta = np.arccos(np.dot(B_vol_iter,PTensor_iter)/(np.linalg.norm(B_vol_iter)*np.linalg.norm(PTensor_iter)))
#   tensor_product = np.array([[B_vol_iter[0]**2, B_vol_iter[0]*B_vol_iter[1], B_vol_iter[0]*B_vol_iter[2]],
#                              [B_vol_iter[0]*B_vol_iter[1], B_vol_iter[1]**2, B_vol_iter[1]*B_vol_iter[2]],
#                              [B_vol_iter[0]*B_vol_iter[2], B_vol_iter[1]*B_vol_iter[2], B_vol_iter[2]**2]])
#   PTensor_iter_rotated = np.dot(I, B_vol_iter)*np.cos(theta) + np.sin(theta) * numpy.cross(B_vol_iter, PTensor_iter) + (1-np.cos(theta))*np.dot(tensor_product, B_vol_iter)
#   PTensor_rotated.append(PTensor_iter_rotated)
#
#
#
#PTensor_rotated = np.array(PTensor_rotated)
#R = rotation_matrix(
#PTensor_rotated = [np.dot(np.dot(rotation_matrix(B_vol[i], np.arccos(B_vol[i][2]/np.linalg.norm(B_vol))),PTensor[i]), rotation_matrix(B_vol[i], np.arccos(B_vol[i][2]/np.linalg.norm(B_vol))).transpose()) for i in xrange(len(B_vol))]

#def rotate( vectorToRotate, vector ):
#   vectorToRotate = v1
#   vector = v2
#   angle = np.arccos(v1.dot(v2)
#   # Get rotation matrix

def rotateTensorToVector( Tensor, vector ):
   '''
      Rotates a tensor with a rotation matrix that would align vector with the z-axis (E.g. moves Tensor to a coordinate system where z axis points in the same direction as vector
      :param Tensor          Tensor to be rotated
      :param vector          Vector for creating the rotation matrix
      :returns rotated tensor
   '''
   vector_u = np.cross(vector, np.array([0,0,1]))
   vector_u = vector_u / np.linalg.norm(vector_u)
   angle = np.arccos( vector.dot(np.array([0,0,1])) / np.linalg.norm(vector) )
   # A unit vector version of the given vector
   R = rotation_matrix( vector_u, angle )
   # Rotate Tensor
   Tensor_rotated = R.dot(Tensor).dot(R.transpose())
   return Tensor_rotated

def rotateVectorToVector( vector1, vector2 ):
   vector_u = np.cross(vector2, np.array([0,0,1]))
   vector_u = vector_u / np.linalg.norm(vector_u)
   angle = np.arccos( vector2.dot(np.array([0,0,1])) / np.linalg.norm(vector2) )
   # A unit vector version of the given vector
   R = rotation_matrix( vector_u, angle )
   # Rotate vector
   vector_rotated = R.dot(vector1)
   return vector_rotated

def rotation_matrix(vector, angle):
   ''' Creates a rotation matrix that rotates into a given vector by a given angle
       :param vector        Some unit vector
       :param angle         Some angle
       :returns a rotation matrix
   '''
   v = vector
   t = angle
   m = np.array([[np.cos(t)+v[0]**2*(1-np.cos(t)), v[0]*v[1]*(1-np.cos(t))-v[2]*np.sin(t), v[0]*v[2]*(1-np.cos(t))+v[1]*np.sin(t)],
                 [v[0]*v[1]*(1-np.cos(t))+v[2]*np.sin(t), np.cos(t)+v[1]**2*(1-np.cos(t)), v[1]*v[2]*(1-np.cos(t))-v[0]*np.sin(t)],
                 [v[0]*v[2]*(1-np.cos(t))-v[1]*np.sin(t), v[2]*v[1]*(1-np.cos(t))+v[0]*np.sin(t), np.cos(t)+v[2]**2*(1-np.cos(t))]])
   return m
