#cython: language_level=3
import numpy as np
#from libc.math cimport csqrtf
#from libc.math cimport cos, sin, sqrt, fmin, fabs, csqrt, cabs
from math import cos, sin, sqrt, sqrt
from cmath import sqrt as csqrt
cimport numpy as np
cimport cython
np.import_array()


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def sdf_vectorized(double t, np.ndarray x,
                   np.ndarray phi):
    cdef double ratio = 1. # length ratio of physical model and numerical model
    cdef double Db = 0.1143 # diameter of PVC trunk
    cdef double Dr = 0.0286 # diameter of PVC root
    cdef np.ndarray [np.float64_t, ndim=1] L = np.array([1.83/ratio,1.83/ratio,1.83/ratio])
    cdef np.ndarray[np.float64_t, ndim=1] theta=np.array(
        [0.0, 0.7853981633974483, 1.5707963267948966, 2.356194490192345, 3.141592653589793, 3.9269908169872414, 4.71238898038469]) 
#    cdef np.ndarray[np.float64_t, ndim=1] X_r=np.array(
#        [0.1595238095238095, 0.21666666666666667, 0.2714285714285714, 0.3285714285714285, 0.38571428571428573, 0.44285714285714284, 0.5])
    cdef np.ndarray[np.float64_t, ndim=1] X_r=np.array(
        [0.67/ratio, 0.91/ratio, 1.14/ratio, 1.38/ratio, 1.62/ratio, 1.86/ratio, 2.1/ratio])
#    cdef np.ndarray[np.float64_t, ndim=1] H_r=np.array(
#        [0.13809523809523808, 0.16904761904761903, 0.19999999999999998, 0.23095238095238094, 0.25952380952380955, 0.29047619047619044, 0.32142857142857145])
    cdef np.ndarray[np.float64_t, ndim=1] H_r=np.array(
        [0.58/ratio, 0.71/ratio, 0.84/ratio, 0.97/ratio, 1.09/ratio, 1.22/ratio, 1.35/ratio])
    cdef int n = phi.shape[0]
    cdef int nr = theta.shape[0]
    cdef int ti, i, ri, gri
    cdef double R,X_0,Y_0,Z_0,A,B,dir_x, dir_y
    cdef double G
#    cdef np.ndarray[np.float64_t, ndim=1] tx = np.array([L[0], L[0], L[0]/2.0,  0.0     , 0.0])
#    cdef np.ndarray[np.float64_t, ndim=1] ty = np.array([0.0 , L[1], L[1]/2.0 , 0.0     , L[1]])
#    cdef np.ndarray[np.float64_t, ndim=1] tx = np.array([0.46/ratio, 1.37/ratio, 0.46/ratio])
#    cdef np.ndarray[np.float64_t, ndim=1] ty = np.array([0.46/ratio, 1.15/ratio, 1.83/ratio])
    cdef np.ndarray[np.float64_t, ndim=1] tx = np.array([0.46/ratio+1.5, 1.37/ratio+1.5, 0.46/ratio+1.5])
    cdef np.ndarray[np.float64_t, ndim=1] ty = np.array([0.46/ratio, 1.15/ratio, 1.83/ratio])
    cdef np.ndarray[np.complex128_t, ndim=1] G_root = np.zeros((3,), dtype=np.complex128)
    cdef double complex num, denom
    for i in range(n):
        phi[i] = 1.0e16
        for ti in range(3):
#        for ti in range(1):
#            ti = 2
            X_0 = x[i,0] - tx[ti]
            Y_0 = x[i,1] - ty[ti]
            Z_0 = x[i,2]
            phi[i] = min(phi[i], sqrt(X_0**2 + Y_0**2) - Db/ratio)
            for ri in range(nr):
                dir_x  = cos(theta[ri]+(ti+1)*2.0943951023931953)#perturb by 2pi/3 so the middle tree has the referene orientation
                dir_y  = sin(theta[ri]+(ti+1)*2.0943951023931953)
                A = -H_r[ri]/X_r[ri]**2
                B = H_r[ri]
                num = 2**(1/3)*3**(2/3)*(3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3)*(1 + csqrt(3)*1.0j)**2 + 12*2**(1/3)*(-2*A*B + 2*A*Z_0 - dir_x**2 - dir_y**2))
                denom = 36*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)*(1 + csqrt(3)*1.0j)
                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
                    G_root[0] = num/denom
                else:
                    G_root[0] = num
                num = 2**(1/3)*3**(2/3)*(3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3)*(1 - csqrt(3)*1.0j)**2 + 12*2**(1/3)*(-2*A*B + 2*A*Z_0 - dir_x**2 - dir_y**2))
                denom = 36*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)*(1 - csqrt(3)*1.0j)
                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
                    G_root[1] = num/denom
                else:
                    G_root[1] = num
                num = 2**(1/3)*3**(2/3)*(-3**(2/3)*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(2/3) + 3*2**(1/3)*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2))
                denom = 18*A**2*((csqrt(3)*A**2*csqrt((27*A**2*(X_0*dir_x + Y_0*dir_y)**2 + 2*(2*A*B - 2*A*Z_0 + dir_x**2 + dir_y**2)**3)/A**6) - 9*X_0*dir_x - 9*Y_0*dir_y)/A**2)**(1/3)
                if abs(denom) > abs(num)*1.0e-16 + 1.0e-16:
                    G_root[2] = num/denom
                else:
                    G_root[2] = num
                for gri in range(3):
                    R = G_root[gri].real
                    if R < -X_r[ri]:
                        R = -X_r[ri]
                    if R > X_r[ri]:
                        R = X_r[ri]
                    G = sqrt((R*dir_x - X_0)**2 + (R*dir_y - Y_0)**2 + (A*R**2 + B - Z_0)**2) - Dr/ratio;
                    phi[i] = min(phi[i], G)
                

#stl_file = "S01T01_mesh.stl"
stl_file = "test2.stl"
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def sdf_vectorized_stl(double t, np.ndarray x, np.ndarray phi):
    cdef int i,j,k,ntri
    with open(stl_file,'r') as f:
        lines = f.readlines()
    ntri = (len(lines)-2)//7
    cdef np.ndarray [np.float64_t, ndim=2] normals = np.zeros((ntri,3),np.float64)
    cdef np.ndarray [np.float64_t, ndim=3] vertices = np.zeros((ntri,3,3),np.float64)
    cdef np.ndarray [np.int32_t, ndim=2] triangles = np.zeros((ntri,3),np.int32)
    nodes_dict={}
    node_tris_dict={}
    #went to a bit of trouble to get triangle neighbors so we can compute an "edge normal"
    #the edge normal is the average of the face normals attached to the edge
    cdef np.ndarray [np.int32_t, ndim=2] triangle_neighbors = -np.ones((ntri,3),np.int32)
    cdef np.ndarray [np.int32_t, ndim=2] triangle_edges = np.zeros((ntri,3),np.int32)
    edges_dict={}
    edge_tris_dict={}
    nN=0
    eN=0
    for i in range(ntri):
        ns = lines[i*7+1].split()[2:]
        for j in range(3):
            normals[i,j] = -float(ns[j])
        for j in range(3):
            vertex = lines[i*7+3+j].split()[1:]
            for k in range(3):
                vertices[i,j,k] = float(vertex[k])
            nt = tuple(vertices[i,j])
            if nt not in nodes_dict:
                nodes_dict[nt] = nN
                node_tris_dict[nN] = []
                nN +=1
            triangles[i,j] = nodes_dict[nt]
            node_tris_dict[triangles[i,j]].append(i)
        for j in range(3):
            nL = triangles[i,j]
            nR = triangles[i,(j+1)%3]
            if nL < nR:
                et = (nL, nR)
            else:
                et = (nR, nL)
            if et not in edges_dict:
                edges_dict[et] = eN
                edge_tris_dict[eN] = []
                eN+=1
            triangle_edges[i,(j+2)%3] = edges_dict[et] #local edge number is local number of opposite node
            edge_tris_dict[edges_dict[et]].append(i)
    cdef np.ndarray [np.float64_t, ndim=2] nodes = np.zeros((nN,3),np.float64)
    for i in range(ntri):
        for j in range(3):
            for k in range(3):
                nodes[triangles[i,j],k] = vertices[i,j,k]
            for k in edge_tris_dict[triangle_edges[i,j]]:
                if i != k:
                    if len(edge_tris_dict[triangle_edges[i,j]]) > 2:
                        print("non-manifold edge", edge_tris_dict[triangle_edges[i,j]], i,k,j)
                        for T in edge_tris_dict[triangle_edges[i,j]]:
                            print(triangles[T])
                    triangle_neighbors[i,j] = k
    from scipy.spatial import KDTree
    tree = KDTree(nodes)
    cdef int nphi = phi.shape[0]
    cdef np.ndarray [np.float64_t, ndim=1] x0 = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] x1 = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] x2 = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] normal = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] edge_normal = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] normal_average = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] normal_average_edge = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] t10 = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] t20 = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] v = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] veMin = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] eMin = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] e = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=1] xHat = np.zeros((3,),np.float64)
    cdef np.ndarray [np.float64_t, ndim=2] A = np.zeros((3,3),np.float64)
    cdef double s, h, hAvg, r, he, heMin, hfMin, sf, se
    for i in range(nphi):
        #get closest node in stl
        dn,nN = tree.query(x[i])
        phi[i] = dn #maximum possible distance to triangulated surface
        s = 1.0
        #fix the sign of phi after we've finished checking for closer entities
        #x[i] may be closer to points on edges or triangles connected to node nN
        #loop over triangles connected to the node
        #accumulate the arithmetic average of the normals of the attached faces
        normal_average[:]=0.0
        heMin=dn #closest edge, if closer than nN
        hfMin=dn #closest face, if closer than nN and any closer edge
        for j in node_tris_dict[nN]:
            normal = normals[j]/np.linalg.norm(normals[j])#ensure it's a unit normal
            normal_average += normal
            x0 = nodes[triangles[j,0]]
            x1 = nodes[triangles[j,1]]
            x2 = nodes[triangles[j,2]]
            v = x[i] - x0
            h = np.dot(v,normal) #distance to plane containing triangle (not necessarily to triangle)
            #x can be closer to a point on the edge of the triangle than to the node
            for k in range(3):
                e = nodes[triangles[j,(k+1)%3]] - nodes[triangles[j,k]] 
                if np.linalg.norm(e) > 0.0:#make sure edge is not degenerate
                    r = - (((nodes[triangles[j,k]]-x[i])*e).sum())/((e*e).sum()) #minimize on line containing edge
                    if 0.0 < r < 1.0: #minimizer is on edge interior
                        ve = x[i] - (nodes[triangles[j,k]]+r*e)#displacement from minimizer
                        he = np.linalg.norm(ve)
                        if he < heMin:
                            heMin = he
                            if he < phi[i]:
                                phi[i] = heMin
                                edge_normal[:] = normal
                                #average with neighbor across edge, if it has one
                                j_neighbor = triangle_neighbors[j,(k+2)%3] 
                                if j_neighbor != -1:
                                    edge_normal += normals[j_neighbor]/np.linalg.norm(normals[j_neighbor])
                                    edge_normal *= 0.5
                                if np.dot(edge_normal, ve) < 0.0:
                                    se = -1.0
                                else:
                                    se = 1.0
            #x[i] can be closer to points inside the face than the node or edge
            t10 = x1-x0
            t20 = x2-x0
            A[:,0] = t10
            A[:,1] = t20
            A[:,2] = normal
            #only check non-degenerate triangles
            if abs(np.linalg.det(A)) > 1.0e-16:
                xHat =  np.linalg.solve(A,v)
                #check that projection into plane is actually in the triangle
                if 0.0 <= xHat[0] <= 1.0 and 0.0 <= xHat[1] <= 1.0 and xHat[0] + xHat[1] <= 1.0:
                    #the zHat should be the same as the projection to plane
                    #assert abs(xHat[2]-h) < 1.0e-16, print(h,xHat[2])
                    #check mapping back into space
                    #assert np.linalg.norm(x0 + np.matmul(A,xHat) - x[i]) < 1.0e-16, print("node and re-projected node", x[i], x0+np.matmul(A,xHat))
                    #assert np.linalg.norm(x0 + xHat[0]*t10 + xHat[1]*t20 + h*normal - x[i]) < 1.0e-16, print("node and re-projected node 2",x[i], x0+xHat[0]*t10+xHat[1]*t20+h*normal)
                    #if the point projected into the plane is in this triangle,
                    #then it must be closer to the triangle than the node or edge
                    #assert abs(h) <= heMin, print(h, heMin)
                    #assert abs(h) <= dn, print(h,dn)
                    if abs(h) < hfMin:
                        hfMin = abs(h)
                        if hfMin < phi[i]:
                            phi[i] = hfMin
                            if h < 0.0:
                                sf = -1.0
                            else:
                                sf = 1.0
        if hfMin < heMin and hfMin < dn:
            s = sf
        elif heMin < dn:
            s = se
        else:
            normal_average /= len(node_tris_dict[nN])
            if np.dot(normal_average, x[i] - nodes[nN]) < 0.0:
                s = -1.0
        phi[i] *= s
