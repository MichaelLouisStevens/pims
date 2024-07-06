import numpy as np
import cmath
cm_sqrt = np.vectorize(cmath.sqrt)

def CircleIntersections(phi, c_offset, lim_aperture, m_offset, mod_aperture, b_offset=0, baf_aperture=1000):

    # if(nargin < 5):
    #     print('insufficient # of input arguments')
    # elif(nargin < 7):
    #     b_offset = 0
    #     baf_aperture = 1000

    c = np.array([c_offset * np.cos(phi), c_offset * np.sin(phi), lim_aperture])
    m = np.array([m_offset * np.cos(phi), m_offset * np.sin(phi), mod_aperture])
    b = np.array([b_offset * np.cos(phi), b_offset * np.sin(phi), baf_aperture])

    n = 3000
    # Need to be divisible by 3
    t = np.linspace(0, 2*np.pi, n+1) 
    c_r = np.sqrt(np.cos(2*t)*(c[0]**2 - c[1]**2) - c[0]**2 + 2*c[0]*c[1]*np.sin(2*t) - c[1]**2 + 2*c[2]**2)/np.sqrt(2) + c[0]*np.cos(t) + c[1]*np.sin(t)
    m_r = cm_sqrt(np.cos(2*t)*(m[0]**2 - m[1]**2) - m[0]**2 + 2*m[0]*m[1]*np.sin(2*t) - m[1]**2 + 2*m[2]**2)/np.sqrt(2) + m[0]*np.cos(t) + m[1]*np.sin(t)
    b_r = cm_sqrt(np.cos(2*t)*(b[0]**2 - b[1]**2) - b[0]**2 + 2*b[0]*b[1]*np.sin(2*t) - b[1]**2 + 2*b[2]**2)/np.sqrt(2) + b[0]*np.cos(t) + b[1]*np.sin(t)

    # converting the complex numbers to absolute
    m_r = np.abs(m_r)
    b_r = np.abs(b_r)

    # clipped region
    min_dis = np.min(np.array([c_r, m_r, b_r]), axis=0)
    # print(np.array([c_r, m_r, b_r]).shape)
    min_ind = np.argmin(np.array([c_r, m_r, b_r]), axis=0)

    # indices of intersection points
    intersect_ind = np.where((min_ind - np.roll(min_ind,-1)) != 0)[0]

    angles_table = np.zeros((len(intersect_ind), 3))
    angles_table = np.append(angles_table, np.array([[1, t[0], min_ind[0]+1],
                                                    [n/3, t[n//3 - 1], min_ind[n//3 - 1]+1],
                                                    [2*n//3, t[2*n//3 - 1], min_ind[2*n//3 - 1]+1],
                                                    [n+1, t[-1], min_ind[-1]+1]]), axis=0)

    # angles_table = np.array([[np.zeros((len(intersect_ind), 3))],
    #                          [0, t[0], min_ind[0]],
    #                          [n/3 - 1, t[n//3 - 1], min_ind[n//3 - 1]],
    #                          [2*n//3 - 1, t[2*n//3 - 1], min_ind[2*n//3 - 1]],
    #                          [n, t[-1], min_ind[-1]]])

    for i in range(len(intersect_ind)):
        angles_table[i,:] = [intersect_ind[i]+1, t[intersect_ind[i]], min_ind[intersect_ind[i]]+1]

    # angles_table = sortrows(angles_table)
    # sorting the rows based on first column
    colsort_ind = np.argsort(angles_table[:,0])
    angles_table = angles_table[colsort_ind,:]

    intersect_dis = min_dis[angles_table[:,0].astype('int')-1]

    return angles_table, min_dis