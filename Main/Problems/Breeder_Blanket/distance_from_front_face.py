


def return_distance_to_surface(a,b,c,d,x,y,z):
  distance = (a * x + b * y + c * z + d ) / ((a*a+b*b+c*c)**0.5)
  print(distance)
  return distance

return_distance_to_surface(a = 0.9908270,
                           b = 0.1158111,
                           c = -0.0696396,
                           d = 1216.1005127,
                           x = 1,
                           y = 1,
                           z = 1)


16.08e6*exp(-11.3*((0.9908270 * x[0] + 0.1158111 * x[1] -0.0696396 * x[2] + 1216.1005127 ) / ((0.9908270*0.9908270+0.1158111*0.1158111+-0.0696396*-0.0696396)**0.5)))
16.08e6*exp(-11.3*(r))

9.96e6*exp(-5.23*((0.9908270 * x[0] + 0.1158111 * x[1] -0.0696396 * x[2] + 1216.1005127 ) / ((0.9908270*0.9908270+0.1158111*0.1158111+-0.0696396*-0.0696396)**0.5)))+18.27e6*exp(-20.56*((0.9908270 * x[0] + 0.1158111 * x[1] -0.0696396 * x[2] + 1216.1005127 ) / ((0.9908270*0.9908270+0.1158111*0.1158111+-0.0696396*-0.0696396)**0.5)))
9.96e6*exp(-5.23*(r))+18.27e6*exp(-20.56*(r))

