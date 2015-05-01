# mdrod

The code in this repository reflects the work done for my senior thesis in applied mathematics at Yale University on generating packings of rod-shaped particles. mdrod102-105 simulate rods of fixed length moving around in a box with no dampening. mdrod120-122 implement the molecular dynamics method to generate packings of the rod-shaped particles. What each program does is explained below.

mdrod102.m
    This program simulates rods moving around in a box where the force on each pair of particles that overlap is calculated as the sum of all the point-to-line or point-to-point distances that is less than 2*R. This produces energy errors when the particles hit side to side or when there are two possible overlaps - the projection and the end-to-end distance both are less than 2*R.

mdrod103.m
    This program uses the same force detector law as mdrod102.m, but here we try to implement a looping method around the center lines of the rods. The ends of the rods are defined by an odd number whereas the line is defined by an even number. So odd-odd interactions are end-to-end distances and odd-even are projections. Even-even interactions are not calculated. Although the number of calculations are the same as mdrod102.m, this program is considerably slower.

mdrod104.m
    This program simulates rods moving around in a box where the force detector only takes the shortest distance calculated from all eight possible distances when it's less than 2*R. This reduces the energy errors by a factor of 10.

mdrod105.m
'This program extends mdrod104.m, but two different size rods can be used. The force detector law is the same as mdrod104.m.
