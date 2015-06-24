from sympy import *



def solve_face_value(integration_limits, has_mdv, has_pdv):
    il=integration_limits
    ans_a=solve([\
                 integrate(rho_hat,(v,il[0][0],il[0][1])) - m3f * (il[0][1] - il[0][0]),\
                 integrate(rho_hat,(v,il[1][0],il[1][1])) - m2f * (il[1][1] - il[1][0]),\
                 integrate(rho_hat,(v,il[2][0],il[2][1])) - m1f * (il[2][1] - il[2][0]),\
                 integrate(rho_hat,(v,il[3][0],il[3][1])) - cf *  (il[3][1] - il[3][0]),\
                 integrate(rho_hat,(v,il[4][0],il[4][1])) - p1f * (il[4][1] - il[4][0]),\
                 integrate(rho_hat,(v,il[5][0],il[5][1])) - p2f * (il[5][1] - il[5][0])],\
                 [a,b,c,d,e,f])
    rho_hat_v=rho_hat.subs([(a,ans_a[a]),(b,ans_a[b]),(c,ans_a[c]),(d,ans_a[d]),(e,ans_a[e]),(f,ans_a[f])])
    rho_hat_ans=rho_hat_v.subs([(v,0)])
    d_rho_hat_ans=simplify(diff(rho_hat_v,v).subs(v,0))
    #print " h6 value "
    #print simplify(rho_hat_ans)
    #print " h5 derivative "
    #print simplify(d_rho_hat_ans)
    if has_mdv and has_pdv:
        for mdv_val in [Rational(1,2), 1, 2]:
            for pdv_val in [Rational(1,2), 1, 2]:
                print( "mdv = ", mdv_val, " pdv = ", pdv_val)
                print( "h6 = ")
                print( simplify(rho_hat_ans.subs([(mdv,mdv_val),(pdv,pdv_val)] )))
                print( "dh5 = ")
                print( simplify(d_rho_hat_ans.subs([(mdv,mdv_val),(pdv,pdv_val)] )))
    elif has_mdv:
        for mdv_val in [Rational(1,2), 1, 2]:
            print( " mdv = ", mdv_val)
            print( " h6 = ")
            print( simplify(rho_hat_ans.subs(mdv,mdv_val)))
            print( " dh5 = ")
            print( simplify(d_rho_hat_ans.subs(mdv,mdv_val)))

    elif has_pdv:
        for pdv_val in [Rational(1,2), 1, 2]:
            print( " pdv = ", pdv_val)
            print( " h6 = ")
            print( simplify(rho_hat_ans.subs(pdv,pdv_val)))
            print( " dh5 = ")
            print( simplify(d_rho_hat_ans.subs(pdv,pdv_val)))

# h6 is the density(-density) function, which when integrated gives us our density values in the cells
# x is the coordinate, which is normalized such that it is x=(v-v_{i-1/2})/dv, where v_{i-1/2} is the velocity at the left face of the center cell i

print( "h6 /dh5 estimates for blockwise AMR grid. cells in current block have dv=1, neighbors are at 1/2 or 2")

a,b,c,d,e,f=symbols('a b c d e f')
v,mdv,pdv=symbols('v mdv pdv')
m3f,m2f, m1f,cf,p1f,p2f=symbols('m3f m2f m1f cf p1f p2f')
rho_hat=a+b*v+c*v**2 + d*v**3 + e*v**4 + f*v**5

print( "-----------------------------------------------------")
print( "First face")
solve_face_value([(-3*mdv,-2*mdv),(-2*mdv,-mdv),(-mdv,0),(0,1),(1,2),(2,3)],True,False)
print( "-----------------------------------------------------")
print( "Second face ")
solve_face_value([(-2*mdv -1 ,-mdv -1),(- mdv-1,-1),(-1,0),(0,1),(1,2),(2,3)],True,False)
print( "-----------------------------------------------------")
print( "Third face ")
solve_face_value([(-mdv -2 ,-2),(-2,-1),(-1,0),(0,1),(1,2),(2,2 + pdv)],True,True)
print( "-----------------------------------------------------")
print( "Fourth face ")
solve_face_value([(-3 ,-2),(-2,-1),(-1,0),(0,1),(1, 1 + pdv ),(1 + pdv ,1 + 2 *pdv)],False,True)
