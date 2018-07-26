Description:
------------
This is an example of the implementation of the control-oriented model for dielectrophoretical force and torque for arbitrarily shaped objects polarized by an arbitrary electric field.
It serves as an electronic supplement to the following publication:
TODO: add citation here

In order to run the example you should have installed the following:
--------------------------------------------------------------------
MATLAB
MATLAB Coder - for faster execution of core components of the model
MATLAB Support for MinGW-w64 C/C++ Compiler
COMSOL Multiphysics with LiveLink (TM) - (LiveLink is necessary for automated basis generation, but it is also possible to run the related cca 20 simulations manually one by one)

The script was tested with the following configurations:
--------------------------------------------------------
1)
MATLAB Version: 9.1.0.441655 (R2016b) with MATLAB Coder and MATLAB Support for MinGW-w64 C/C++ Compiler
Java Version: Java 1.7.0_60-b19 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
COMSOL 5.1 with MATLAB (LiveLink (TM))
Wolfram Mathematica 10.4

Microsoft Windows 10 Home Version 10.0 (Build 17134), 64 bit
Intel(R) Core(TM) i5 CPU M 430 @ 2.27GHz
8GB RAM

2)
MATLAB Version: 9.4.0.813654 (R2018a) with MATLAB Coder and MATLAB Support for MinGW-w64 C/C++ Compiler
Java Version: Java 1.8.0_144-b01 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
COMSOL Multiphysics 5.3a with MATLAB
Wolfram Mathematica 11.2

Microsoft Windows 7 Enterprise  Version 6.1 (Build 7601: Service Pack 1), 64 bit
Intel(R) Core(TM) i5-3550 CPU @ 3.30GHz
8GB RAM

Instructions:
-------------
1. Run the "model_comparision.m" script section by section and follow the included instructions.
   The main parts are:
   - Setting the inputs of the model;
   - The setup script will generate the basis of multipolar moments (as described in the paper) and lookup-tables for potential derivatives and compile the core components of the model for your system;
   Some of these steps might be automatically skipped since the demo already includes some precomputed basis. Generally, this operation might took a significant amount of time (several hours or even units of days on a convetional PC);
   - Computation of the force and torque acting on the object of interest using both the presented control-oriented model and the reference MST simulation;
   - Comparison of the two sets of results;
2. Feel free to inspect the model and try it with various objects of different shapes or placed even above different electrode arrays (There are two different tetris shapes and two different layouts of quadrupolar electrode arrays to try for the beginning.).

Licence and citing:
-------------------
Copyright (C) 2018  Tomas Michalek

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

If you use this piece of software in your research, please cite the following paper
(TODO: fill in the citing info)

Versions:
---------
v1 - initial version containing the example of implementation of a control-oriented EM method for arbitrarily shaped objects

Contact:
--------
Tomas Michalek
PhD student
Czech Technical University in Prague, Faculty of Electrical Engineering
Department of Control Engineering

Email: michato4@fel.cvut.cz
Phone: +420-22435-7471

Address:
Karlovo namesti 13, 121 35, Prague, Czech republic. 

