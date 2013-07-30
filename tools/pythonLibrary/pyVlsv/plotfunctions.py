import numpy as np
from math import pi

functionList = dict()

def list_fit_functions():
   global functionList
   print "Functions:"
   for i in functionList.iteritems():
      print i[0]

class PlotFunction:
   parameters = 1
   def function(self, args, x, y):
      print "Called PlotFunction's function in plotfunctions.py (not allowed)"
   def get_function(self):
      return self.function
   def get_parameters(self):
      return self.parameters

################################################################
#FUNCTION: FITTING A + B*COS(X) + C*SIN(X)
class sincosfit(PlotFunction):
   parameters = 3
   fitLength = 2*pi
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*np.cos(2*pi*x/self.fitLength) + args[2]*np.sin(2*pi*x/self.fitLength)
      return y - value
################################################################

################################################################
#FUNCTION: FITTING A + B*X + C*X^2 + D*X^3
class polynomialthirdfit(PlotFunction):
   parameters = 4
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*x + args[2]*(x**2) + args[3]*(x**3)
      return y - value
################################################################

################################################################
#FUNCTION: FITTING A + B*X + C*X^2
class polynomialsecondfit(PlotFunction):
   parameters = 3
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*x + args[2]*(x**2)
      return y - value
################################################################

################################################################
#FUNCTION: FITTING A + B*X
class polynomialfirstfit(PlotFunction):
   parameters = 2
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*x
      return y - value
################################################################

################################################################
#FUNCTION: FITTING A + B*e^x
class exponentialfirstfit(PlotFunction):
   parameters = 2
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*np.exp(x)
      return y - value
################################################################

################################################################
#FUNCTION: FITTING A + B*e^x + C*e^(-x)
class exponentialsecondfit(PlotFunction):
   parameters = 3
   def function(self, args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*np.exp(x) + args[2]*np.exp((-1)*x)
      return y - value
################################################################

################################################################
#FUNCTION: FITTING EMPTY
class emptyfit(PlotFunction):
   parameters = 1
   def function(self, args, x, y):
      return y
################################################################

################################################################
#LIST FUNCTIONS AND THE NUMBER OF ARGUMENTS THEY HAVE HERE
functionList["sincos"] = sincosfit()
functionList["polynomialthird"] = polynomialthirdfit()
functionList["polynomialsecond"] = polynomialsecondfit()
functionList["polynomialfirst"] = polynomialfirstfit()
functionList["exponentialfirst"] = exponentialfirstfit()
functionList["exponentialsecond"] = exponentialsecondfit()
functionList["emptyfit"] = emptyfit()
################################################################



