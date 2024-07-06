import sys
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
plt.ion()
import CalculateTriantArea

# Generate Lookup Table

# angular resolution of the table
ang_res = 2.5
# theta_max = 44.8716
theta_max = 54.950

theta = np.arange(ang_res, np.floor(theta_max/ang_res)*ang_res + ang_res, ang_res)
phi = np.arange(0, 360-ang_res, ang_res)
# phi = np.array([0, 90, 180, 270])
Area = np.zeros((len(theta) * len(phi), 5))

for i in range(len(theta)):
    for j in range(len(phi)):
        print(theta[i],phi[j])
        AreaA, AreaB, AreaC = CalculateTriantArea.CalculateTriantArea(theta[i], phi[j], 0, 0)
        totalArea = AreaA + AreaB + AreaC
        Area[(i*len(phi)+j), :] = np.array([theta[i], phi[j], AreaA/totalArea, AreaB/totalArea, AreaC/totalArea])

'''

save LUT_0.25_theta_max_55.mat
keyboard
load LUT_0.25_theta_max_55.mat;

%% Plot Solution

figure
scatter3(Area(:,1),Area(:,2),Area(:,3),'.', 'DisplayName', 'Triad A')
hold on
scatter3(Area(:,1),Area(:,2),Area(:,4),'.', 'DisplayName', 'Triad B')
[X,Y] = meshgrid(theta,phi);scatter3(Area(:,1),Area(:,2),Area(:,5),'.', 'DisplayName', 'Triad C')
axis tight
hold off
xlabel('Theta (deg)');
ylabel('Phi (deg)');
zlabel('Area (%)');
legend('show', 'Location','northwest');

figure
[X,Y] = meshgrid(theta,phi);
tri = delaunay(X,Y);
Z = abs(a-Area(2:end,3))+abs(b-Area(2:end,4))+abs(c-Area(2:end,5));
Z_reshaped = reshape(Z,length(theta),length(phi));
 trimesh(tri,X,Y,Z_reshaped)


[val_min, ind_min] = min(Z);
theta_min = Area(ind_min,1);
phi_min = Area(ind_min,2);
'''

# loading the matlab results and overplotting
Area_M = scipy.io.loadmat('../Area_Matlab.mat')['Area']
Area_M = Area_M[1:]

# making the 3D scatter plots
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(Area[:,0], Area[:,1], Area[:,2], '.', s=1, label='Triad A')
ax.scatter(Area[:,0], Area[:,1], Area[:,3], '.', color='red', s=1, label='Triad B')
ax.scatter(Area[:,0], Area[:,1], Area[:,4], '.', color='orange', s=1, label='Triad C')
ax.set_xlabel('Theta (deg)')
ax.set_ylabel('Phi (deg)')
ax.set_zlabel('Area (%)')
ax.set_zlim([0,1])
plt.legend()
# N = -1

# ax = fig.add_subplot()
# ax.plot(Area_M[:N,0], Area_M[:N,2], '.')
# ax.plot(Area_M[:N,0], Area_M[:N,3], '.')
# ax.plot(Area_M[:N,0], Area_M[:N,4], '.')

# ax.plot(Area[:N,0], Area[:N,2], 'x')
# ax.plot(Area[:N,0], Area[:N,3], 'x')
# ax.plot(Area[:N,0], Area[:N,4], 'x')