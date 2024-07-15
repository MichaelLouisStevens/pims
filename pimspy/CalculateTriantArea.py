import numpy as np
import cmath
cm_sqrt = np.vectorize(cmath.sqrt)
from scipy.integrate import quad as integral
import matplotlib.pyplot as plt
plt.ion()

import CircleIntersections as CircInt

def CalculateTriantArea(theta, phi, HV_DC, HV_AC, plot=True):
    # Calculate S-curve Deflection

    deflection = 0.000026851 * theta**3 + 0.0010242 * theta**2 + 0.24561 * theta - 0.029726

    # variance evaluates to NaN in the current setup. Skipping its computation
    # variance = -0.1422036 * np.log(HV_AC/HV_DC) + 0.2692889
    # variance = variance * deflection

    # Cup Geometrc2 - Verifc2 Valid Angles
    phi = phi * np.pi/180
    theta = theta * np.pi/180

    # Cup opening (baffle aperture) radius
    baf_aperture = 200/2
    # Baffle depth
    baf_depth = 38.1
    # Spacing around HV
    lv_space = 7
    # Spacing around LV
    hv_space = 9
    # Modulator Aperature radius
    mod_aperture = 132/2
    # Modulator depth
    # mod_depth = 1.6 + 1 + 2.5 + 1 + 2 * lv_space + 2 * hv_space
    # Modulator depth
    mod_depth = 38.4
    # mod_depth = 1.6+1+2.5+1+2*lv_space+2*hv_space;			% Modulator depth original

    # Limiting Aperature radius
    lim_aperture = 80/2
    # col_plate = 97.8/2;			% Collector plate radius 
    col_plate = 100.8/2
    # col_depth = 9.15;      % FM baseline + shim
    # col_depth = 13.15;      # new2
    col_depth = 11.15;        # new1
    # col_depth = 8.4;
    # col_depth = 8.7;			% Collector depth original

    # hv_cp_distance = 1.6+1+2*lv_space+hv_space+col_depth;
    # baffle_lim_ang = atan((baf_aperture-lim_aperture-deflection)/(baf_depth+mod_depth));     % Baffle Aperture to Limiting Aperture Angle
    # mod_lim_ang = atan((mod_aperture-lim_aperture-deflection)/(mod_depth));                        % Modulator Aperture to Limiting Aperture Angle
    # baffle_cp_center_ang = atan((baf_aperture)/(baf_depth+mod_depth+col_depth));  % Baffle Aperture to center of CP
    # mod_cp_center_ang = atan((mod_aperture)/(mod_depth+col_depth));                     % Modulator Aperture to center of CP
    # full_ang = min([baffle_lim_ang, mod_lim_ang, baffle_cp_center_ang, mod_cp_center_ang]);        % Angles < full_angle --> Whole circle
    # full_ang = 25.9123*pi/180;

    full_ang = 25.79 * np.pi / 180
    # valid_ang = atan((baf_aperture-deflection)/(baf_depth+mod_depth+col_depth))
    valid_ang = 54.9950 * np.pi / 180
    # valid_ang = 44.8716*pi/180
    switchover_offset = (-baf_aperture * mod_depth + baf_depth * mod_aperture\
                        + mod_aperture * mod_depth)/baf_depth
    switchover_angle = np.arctan((baf_aperture - switchover_offset) / (baf_depth + mod_depth))

    # Calculate Area
    c_offset = col_depth * np.tan(theta)
    m_offset = (col_depth + mod_depth) * np.tan(theta) + deflection
    b_offset = (col_depth + mod_depth + baf_depth) * np.tan(theta) + deflection
    c1 = c_offset * np.cos(phi)
    c2 = c_offset * np.sin(phi)
    m1 = m_offset * np.cos(phi)
    m2 = m_offset * np.sin(phi)
    b1 = b_offset * np.cos(phi)
    b2 = b_offset * np.sin(phi)

    def c_fun(psi):
        return (cm_sqrt(np.cos(2 * psi) * (c1**2 - c2**2) - c1**2 + 2*c1*c2 * np.sin(2.*psi) - c2**2 + 2*lim_aperture**2) / np.sqrt(2) + c1 * np.cos(psi) + c2 * np.sin(psi))**2

    def m_fun(psi):
        return (cm_sqrt(np.cos(2 * psi) * (m1**2 - m2**2) - m1**2 + 2*m1*m2 * np.sin(2.*psi) - m2**2 + 2*mod_aperture**2) / np.sqrt(2) + m1 * np.cos(psi) + m2 * np.sin(psi))**2

    def b_fun(psi):
        return (cm_sqrt(np.cos(2 * psi) * (b1**2 - b2**2) - b1**2 + 2*b1*b2 * np.sin(2.*psi) - b2**2 + 2*baf_aperture**2) / np.sqrt(2) + b1 * np.cos(psi) + b2 * np.sin(psi))**2

    t = np.linspace(0, 2 * np.pi, 1000)
    CP = np.array([col_plate * np.sin(t), col_plate * np.cos(t)])
    triant = np.array([[0, col_plate, 0, col_plate*np.cos(2*np.pi/3), 0, col_plate*np.cos(4*np.pi/3)],
                       [0, 0, 0, col_plate*np.sin(2*np.pi/3), 0, col_plate*np.sin(4*np.pi/3)]])
    proj_c = np.array([c1+lim_aperture*np.sin(t), c2+lim_aperture*np.cos(t)])
    proj_m = np.array([m1+mod_aperture*np.sin(t), m2+mod_aperture*np.cos(t)])
    proj_b = np.array([b1+baf_aperture*np.sin(t), b2+baf_aperture*np.cos(t)])
    '''
    plot(CP(1,:),CP(2,:),'r',triant(1,:),triant(2,:), 'r-', proj_c(1,:),proj_c(2,:),'g.-',proj_m(1,:),proj_m(2,:),'k.-',proj_b(1,:),proj_b(2,:),'b.-')
    # plot(CP(1,:),CP(2,:),'r',triant(1,:),triant(2,:), 'r-', proj_c(1,1:600),proj_c(2,1:600),'k',proj_c(1,730:end),proj_c(2,730:end),'k',proj_m(1,600:735),proj_m(2,600:735),'k')
    axis equal
    '''

    # Signal not on all 3 CP (with deflection) - >44.871
    if(theta > valid_ang):
        AreaA = -1
        AreaB = -1
        AreaC = -1

    # Full Circle Projected - <25.913
    elif(theta < full_ang):
        AreaA = integral(c_fun, 0, 2*np.pi/3)[0] / 2
        AreaB = integral(c_fun, 2*np.pi/3, 4*np.pi/3)[0] / 2
        AreaC = integral(c_fun, 4*np.pi/3, 2*np.pi)[0] / 2

    # Projection clipped by Modulator Aperture
    elif(theta < switchover_angle):
        angles_table, __ = CircInt.CircleIntersections(phi, c_offset, lim_aperture, m_offset, mod_aperture, b_offset, baf_aperture)

        AreaA = 0
        AreaB = 0
        AreaC = 0

        # padding an extra row at the end
        angles_table = np.pad(angles_table, (0,1))[:-1]

        for i in range(len(angles_table) - 1):
            if(angles_table[i+1, 2] == 1):
                angles_table[i, 3] = integral(c_fun, angles_table[i, 1], angles_table[i+1, 1])[0] / 2
            else:
                angles_table[i, 3] = integral(m_fun, angles_table[i, 1], angles_table[i+1, 1])[0] / 2
            
            # Triad A
            if(angles_table[i+1, 1] < 2 * np.pi / 3):
                AreaA = AreaA + angles_table[i,3]
            # Triad B
            elif(angles_table[i+1, 1] < 4 * np.pi / 3):
                AreaB = AreaB + angles_table[i, 3]
            # Triad C
            else:
                AreaC = AreaC + angles_table[i, 3]
        
    # Projection clipped by Baffle Aperture and Modulator Aperture
    else:
        angles_table, __ = CircInt.CircleIntersections(phi, c_offset, lim_aperture, m_offset, mod_aperture, b_offset, baf_aperture)
        
        AreaA = 0
        AreaB = 0
        AreaC = 0

        # padding an extra row at the end
        angles_table = np.pad(angles_table, (0,1))[:-1].astype('complex128')
        
        for i in range(len(angles_table)-1):
            if(angles_table[i+1, 2] == 1):
                angles_table[i, 3] = integral(c_fun, angles_table[i, 1], angles_table[i+1, 1], complex_func=True)[0]/2
            elif(angles_table[i+1, 2] == 2):
                angles_table[i, 3] = integral(m_fun, angles_table[i, 1], angles_table[i+1, 1], complex_func=True)[0]/2
            else:
                angles_table[i, 3] = integral(b_fun, angles_table[i, 1], angles_table[i+1, 1], complex_func=True)[0]/2
            
            # Triad A
            if(angles_table[i+1,1] < 2*np.pi/3):
                AreaA = AreaA + angles_table[i,3]
            # Triad B
            elif(angles_table[i+1,1] < 4*np.pi/3):
                AreaB = AreaB + angles_table[i,3]
            # Triad C
            else:
                AreaC = AreaC + angles_table[i,3]

    '''
    if(plot):
        m_offset = (col_depth+mod_depth)*tan(theta)+variance
        b_offset = (col_depth+mod_depth+baf_depth)*tan(theta)+variance
        m1 = m_offset*cos(phi)
        m2 = m_offset*sin(phi)
        b1 = b_offset*cos(phi)
        b2 = b_offset*sin(phi)
        proj_m = [m1+mod_aperture*sin(t);m2+mod_aperture.*cos(t)]
        proj_b = [b1+baf_aperture*sin(t);b2+baf_aperture.*cos(t)]
        hold on
        plot(proj_m(1,:),proj_m(2,:),'k--',proj_b(1,:),proj_b(2,:),'b--')
        % plot(proj_m(1,630:710),proj_m(2,630:710),'k')
        hold off

        sum = AreaA+AreaB+AreaC

        plt.figure()
        axis([-60, 60, -60, 60])
        plt.text(25, 50, f'Area A = {AreaA}', color='black', fontsize = 14)
        plt.text(30, 47, num2str(AreaA/sum*100) '%', color='black', fontsize = 14)
        plt.text(-55, 50, f'Area B = {AreaB}', color='black', fontsize = 14)
        plt.text(-50, 47, num2str(AreaB/sum*100) '%', color='black', fontsize = 14)
        plt.text(25, -47, f'Area C = {AreaC}', color='black', fontsize = 14) # testing 
        plt.text(30, -50, num2str(AreaC/sum*100) '%', color='black', fontsize = 14)
    '''

    return AreaA, AreaB, AreaC
