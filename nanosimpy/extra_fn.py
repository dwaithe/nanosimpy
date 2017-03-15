def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > ((B[1]-A[1]) * (C[0]-A[0]))

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)
from __future__ import division
def dist_meas(R,Ax,Ay,Bx,By,Cx,Cy):
    Q = np.array([Cx,Cy])            # Centre of circle
    r = R                 # Radius of circle
    P1 = np.array([Ax,Ay])      # Start of line segment
    V = np.array([Bx,By]) - P1  # Vector along line segment
    a = V.dot(V)
    b = 2 * V.dot(P1 - Q)
    c = P1.dot(P1) + Q.dot(Q) - 2 * P1.dot(Q) - r**2
    disc = b**2 - 4 * a * c
    #Doesn't intersect
    if disc < 0:
        
        return False#, None
    sqrt_disc = math.sqrt(disc)
    t1 = (-b + sqrt_disc) / (2 * a)
    t2 = (-b - sqrt_disc) / (2 * a)
    if np.sqrt((Ax-Cx)**2 + (Ay-Cy)**2) < R and np.sqrt((Bx-Cx)**2 + (By-Cy)**2) < R:
            return True
    if not (0 <= t1 <= 1 or 0 <= t2 <= 1):
        return False#, None
    t = max(0, min(1, - b / (2 * a)))
    return True#, P1 + t * V
def test_all_intersect( x1,i,vector_list):
    Ax = x1[0,i-1]
    Bx = x1[0,i]
    Ay = x1[1,i-1]
    By = x1[1,i]
    #plot(Ax,Ay,'o')
    dist = []
    dist_arr = []
    an_intersect = False
    tailored_list = []
    
    Axh = Ax - (Ax%200)
    Bxh = Ax - (Ax%200)
   
    tailored_list = vector_list[Axh,Bxh]
    
        
    for vec in tailored_list:
        Cx = narrow_list[vec][0]
        Cy = narrow_list[vec][2]
        Dx = narrow_list[vec][1]
        Dy = narrow_list[vec][3]
        E = narrow_list[vec][4]
        plot([Cx,Dx],[Cy,Dy],'k')
        if intersect([Ax,Ay],[Bx,By],[Cx,Cy],[Dx,Dy]) == True:
            an_intersect = True
            #Calculate the point of intersection.
            denom = ((Ax-Bx)*(Cy-Dy)-(Ay-By)*(Cx-Dx))
            iX = (((Ax*By - Ay*Bx)*(Cx-Dx)-(Ax-Bx)*(Cx*Dy - Cy*Dx))/denom)
            iY = (((Ax*By - Ay*Bx)*(Cy-Dy)-(Ay-By)*(Cx*Dy - Cy*Dx))/denom)
            
            
            dist.append(np.sqrt((iX-Ax)**2 + (iY-By)**2))
            dist_arr.append([Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,x1,iX,iY,i,E])
    
    if an_intersect == True:
        it = np.argmin(dist)
        
        if dist_arr[it][12] == 1:
            print '1'
        if np.random.random() > 1.1 or dist_arr[it][12] == 1:
                return
        
        x1 = calculate_norm_ray(*dist_arr[it])
        test_all_intersect(x1,i,vector_list)
def calculate_norm_ray(Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,x1,iX,iY,i,E):
    normalY = Dx - Cx
    normalX = Cy - Dy
    normalLen = np.sqrt(normalX**2 + normalY**2)
    normalX = normalX/normalLen
    normalY = normalY/normalLen

    rayX = Bx - iX
    rayY = By - iY
    dotProduct = (rayX*normalX)+(rayY*normalY)
    dotNormalX = dotProduct * normalX
    dotNormalY = dotProduct * normalY

    reflectedRayTipX = Bx - (dotNormalX*2)
    reflectedRayTipY = By - (dotNormalY*2)

    #Replaces last value with new coordinate.
    x1[:,i] = [reflectedRayTipX,reflectedRayTipY]
    
    return x1


def brownian(n, dt, delta, coord,vec_list,num_of_mol):
   
    
    full_arr = {}
    for b in range(0,num_of_mol):
        full_arr[b] = np.zeros((2,N))
        full_arr[b][:,0] = coord[:,b]
    for i in range(1, N):
        r = norm.rvs(size=[2,num_of_mol], scale=delta*sqrt(dt))
        for b in range(0,num_of_mol):
            full_arr[b][:,i] += full_arr[b][:,i-1] + r[:,b]
        
       
            test_all_intersect(full_arr[b],i,vector_list)
    return  full_arr