#General

This is a project repository for physically based simulation class.

#Projects

1.Example:       opengl viewer test project
2.BalloonTest:   balloon simulation test project
3.WaterTest:     SPH simulation test project (spherical boundary)
4.WaterCubeTest: SPH simulation test project (cubic boundary)
5.WaterBalloon:  water balloon simulation project (final project)
6.Optimization:  SPH accelertion project

#How to build

1. Download Eigen and place it in the project root directory (ex. PBSProject/Eigen/).
2. Download libigl and place it in the project root directory (ex. PBSProject/libigl/).
3. Compile libigl.
4. Go to a project directory and create a build folder, and use cmake to create make file.
5. Make it!

#How to use

##Key

e: scaling view mode
r: rotation view mode
t: change drawing mode of triangles
p: play or stop
h: switch pomp (air or water)

##WaterBalloon

Changing "useWater" flag in main.cpp, you can switch whether you simulate air or fluids.
