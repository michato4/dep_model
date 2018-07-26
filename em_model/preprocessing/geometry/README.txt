The object should be of similar size as is the "tetris_s" object (It should fit into a ball having radius of 95e-6m and it should not be much smaller). Otherwise a modification of the Comsol model "extraction_of_multipoles.mph" will be necessary.
Furthermore, to avoid numerical issues during FEM processing, round at least a little bit all the edges of the object.

Similarly, the size of the electrodes should be similar to the two examples contained in this folder.
Furhtermore, the arrangement of the electrodes should be similar (symmetric quadrupolar arrangement), since otherwise further modifications of the Comsol model "exact_boundary_conditions.mph" and the "setup.m" script will be necessary.

These requirements were posed to keep the example scripts as simple as possible, but they do not represent limitations of the method itself.

The objects were created in Inventor, but any program supporting export into .stl can be used.
The electrodes were drawn in AutoCAD, but any program supporting export into .dxf can be used.