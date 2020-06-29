import matplotlib.pyplot as plt

fig, ax = plt.subplots()
fig.set_size_inches(10, 10)


#original detections
circle = plt.Circle((8.8725, -10.5872), 18.56/60/2, color='r', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.15, -10.5872), 18.56/60/2 , color='r', fill=False)
ax.add_artist(circle)

# change default range so that new circles will work
ax.set_aspect('equal')

ax.set_xlim((8.71, 9.31))
ax.set_ylim((-10.29, -10.89))
plt.xlabel(r"Right Acension $^\circ$", fontsize=20)
plt.ylabel(r"Declination $^\circ$", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='minor', labelsize=20)
# some data
fig.savefig('original_detections.png', dpi=100, bbox_inches = 'tight')

#P2C follow up
circle = plt.Circle((9.0113, -10.5872), 18.56/60/2, color='r', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.0113, -10.7259), 18.56/60/2 , color='r', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.0113, -10.4485), 18.56/60/2 , color='r', fill=False)
ax.add_artist(circle)

fig.savefig('P2C_detections.png', dpi=100, bbox_inches = 'tight')

#P1 follow up
circle = plt.Circle((9.06, -10.59),  2.53/60/2, color='b', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.039, -10.59),  2.53/60/2, color='b', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.081, -10.59),  2.53/60/2, color='b', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.06, -10.611),  2.53/60/2, color='b', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.06, -10.569),  2.53/60/2, color='b', fill=False)
ax.add_artist(circle)

fig.savefig('P1_detections.png', dpi=100, bbox_inches = 'tight')

#P2E follow up
circle = plt.Circle((9.059, -10.587),  1.26/60/2, color='g', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.059, -10.577),  1.26/60/2, color='g', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.059, -10.597),  1.26/60/2, color='g', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.069, -10.587),  1.26/60/2, color='g', fill=False)
ax.add_artist(circle)
circle = plt.Circle((9.049, -10.587),  1.26/60/2, color='g', fill=False)
ax.add_artist(circle)

fig.savefig('P2E_detections.png', dpi=100, bbox_inches = 'tight')
