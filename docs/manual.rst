.. _Fener_manual:


***************
Fener Manual
***************

.. _introduction:

Introduction
=============================

Fener evaluates the thermal, daylighting and glare performance of 
complex fenestration systems in office spaces.

Fener uses the Radiance-based three-phase method and bi-directional 
scattering distribution functions (BSDF) to make irradiance and 
illuminance calculations. Glare analysis is performed based on the 
Simplified DGP and Enhanced Simplified DGP methods (the latter is only 
implemented at developer level). Thermal calculations are carried out 
according to the Kuhn2011 model for the heat transfer through the 
window, which additionally requires a U-value and Directional Solar Heat
Coefficients (DSHGC) as the thermal characterization of a fenestration 
system. Transfer functions and the energy balance method make it 
possible to calculate indoor air temperature and energy loads. Some 
auxiliary programs are also available to calculate BSDF and DSHGC for 
simplified composed of a shading device and a glazing unit. Fener is 
also able to generate an input data file for EnergyPlus and run it, 
using the results of the daylighting calculation. 

.. _examples:

Examples
=============================

For a completely new case study, where also the weather data file is 
used for the first time, a full Fener command would look like::

  > python3 ./src/fener.py {directory_path`/config.fnr -c 20 -meteo -geo -shoeBox -daylight -grid -lightSch -therm -glareSimpl

This case study consists of 20-core simulation from a shoe-box room 
geometry. The program generates the meteorological files, the Radiance 
geometry (including a grid of horizontal photosensors) and the 
three-phase method matrices. The program runs daylight, glare and 
thermal simulations. A schedule of artificial light operation is 
generated and used in the thermal simulation. The program generates an 
input data file (idf) for EnergyPlus and runs it.

Once the meteorological and geometry files are created, running 
variations of the same case study that share these files can be realized 
by running::

  > python3 ./src/fener.py {directory_path`/config.fnr -daylight -grid -lightSch -therm -glareSimpl

On the other hand, if the goal is to perform just a daylight simulation 
for a grid of horizontal photosensors (no thermal or glare 
calculations), the command line is the one below. In this case, 
meteorological files, room geometry and three-phase method matrices have 
already been created from previous simulations::

  > python3 ./src/fener.py {directory_path`/config.fnr -daylight -grid

.. _program_options:

Program Options
=============================

In the program commands above, :file:`fener.py` is the main routine of the program, which contains calls to all other routines, and :file:`config.fnr` is the input data file. The rest are options of the program described below:

Option :file:`-c` sets the number of cores (processes) to accelerate computation on a shared memory machine. This option is highly recommended to be used when using Radiance (e.g. generation of three-phase method matrices (option :file:`-geo`) or enhanced-DGP glare calculations (option :file:`-glare`).

Option :file:`-meteo` generates required files from the meteorological file. This option must be used when a new meteorological data file is used. Generated files are stored in the same folder as the meteorological file and are available for other simulation with the same meteorological data.

Option :file:`-geo` generates the Radiance geometry and three-phase method matrices of the case study from the program inputs. This option must be used when a new geometry (or a change in the geometry) is specified. The files generated are stored in the folder :file:`workDir/` and are available for other simulations with the same geometry.

Option :file:`-shoeBox` builds a shoe-box geometry (rotatable rectangular space). This option creates a file :file:`surfaces.dat`, which has the multidimensional inputs of opaque surfaces.

Option :file:`-grid` builds a horizontal grid of illuminance sensors. This option must be used combined with the options :file:`-geo` and :file:`-daylight` or with :file:`-geo` and :file:`-df`. 

Option :file:`-display` visualizes the geometry using :file:`objline`, radiance quick preview program. Simulation is paused during displaying. This option must be used combined with the option :file:`-geo`.

Option :file:`-daylight` runs a daylighting simulation.
Without the :file:`-grid` option, a Radiance description of horizontal illuminance sensors must be provided and referred in :file:`config` file.

Option :file:`-lightSch` calculates a schedule of artificial lighting operation from the daylight simulation.

Option :file:`-glareSimpl` runs a glare simulation based on the Simplified-DGP method. The Simplified DGP is proportional to the vertical illuminance at predefined sensor positions. 
A Radiance description of vertical illuminance sensors must be provided and referred in the :file:`config.fnr` file.

Option :file:`-glare` runs a glare simulation based on the Enhanced Simplified DGP method. The glare model requires a Radiance simulation every timestep, which implies a added computational cost to the simulations. A Radiance geometry and material definition of the fenestration system must be provided and referred in the :file:`config.fnr` file.

Option :file:`-therm`runs a thermal simulation. The transient heat transfer through opaque elements (walls, ceiling) is solved based on transfer functions. The heat transfer through fenestration systems is solved based on the Black-Box model (Kuhn et al 2011).
Angular-dependent g-values or Directional Solar Heat Coefficients (DSHGC) of the fenestration system must be provided and referred in the :file:`config.fnr` file.

Option :file:`-thermComf`, based on the computer model of the standard ASHRAE 55-2010 for thermal comfort. This module calculates a seven-point thermal sensation scale (PMV) and a percentage of dissatisfied occupants (PPD) (ANSI/ASHRAE Standard 55-2010).
The option :file:`-thermComf` can only be run together with :file:`-therm`.

Option :file:`-df` runs a daylighting factor simulation. The sky distribution for this calculation corresponds to a standard CIE overcast day. 

Option :file:`-da` performs a daylight autonomy calculation. This option must be used combined with the option :file:`-daylight`.

Option :file:`-klems` generates a BSDF dataset and angular-dependent layer absortivities of a fenestration system from the BSDF of the individual layers based on the Klems' method. For glazing and perfect diffusive layers (shades), \textbf{Fener` can calculate the BSDF internally. BSDF files are saved according to the paths indicated in the :file:`config.fnr` file. 

Option :file:`-outside` includes an outdoor scene in Radiance format in the Daylight matrix calculation. A file containing the outdoor scene must be specified in the :file:`config.fnr` file. 

Option :file:`-schCntrl` controls a fenestration system based on a predefined schedule of window construction. The schedule must be specified in the :file:`config.fnr` file.

Option :file:`-mtxCntrl` controls a fenestration system based on a matrix of conditions that depend on internal variables (e.g. indoor temperature, solar radiation on the facade, etc.).

.. _input-data-file:

Input Data File
=============================

The :file:`config` file is the input data file of Fener. Config files are located in the directory :file:`config.fnr` in the main folder of the program. It is composed of two sections: :file:`PATHS` and :file:`VARIABLES`. It may include comments, prefixed by the character #. Unless otherwise indicated, units are always IS (e.g. m, K, W, etc.). The inputs required by each program option are indicated below.

:file:`PATHS`
This section of the configuration file tells the program which directories to search for files in which multidimensional inputs are defined. In \textbf{Fener}, multidimensional inputs are always defined in text files. For example, the construction of a wall is a multidimensional input variable that contains the material indexes for each layer of the wall. Multidimensional inputs are the following:

:file:`meteo` 
Meteorological information: air temperature, direct normal solar irradiation, diffuse horizontal solar irradiation, incoming longwave radiation and wind speed.

:file:`input`
Folder for simulation inputs.

:file:`output`
Folder for simulation outputs.

:file:`workDir`
Folder for simulation working files.

:file:`matOpaq`
Opaque material database.

:file:`constOpaq`
Construction of opaque areas.

:file:`surf`
Geometry of the enclosure.

:file:`win`
Geometry of glazing areas.

:file:`frame`
Geometry of frames.

:file:`bsdfSys_0`
BSDF dataset of the fenestration system (as many as window constructions). This information can be generated with option :file:`-klems` or externally with the program WINDOW.

:file:`illuPts`
Radiance description of horizontal illuminance sensors. Not necessary if :file:`-grid`

:file:`lightCntrl`
Input parameters for the daylight control of artificial lighting.

:file:`lightSch`
Operation schedule of artificial lighting.

:file:`calorim_0`
Calorimetric dataset of the fenestration system. This information can be generated with option :file:`-calorim`

:file:`infSch`
Operation schedule of infiltration.

:file:`occSch`
Operation schedule of occupation.

:file:`equipSch`
Operation schedule of equipment.

:file:`heatSpSch`
Operation schedule of heating thermal setpoints.

:file:`coolSpSch`
Operation schedule of cooling thermal setpoints.

:file:`winRad_0`
Radiance geometry and material of the fenestration system.

:file:`illuVertPts`
Radiance description of vertical illuminance sensors.

:file:`matGlz`
Glazing material database.

:file:`matGas`
Window gas material database.
 
:file:`matBSDF`
BSDF layer material database (thermal description).
 
:file:`constWin`
Construction of glazing areas.

:file:`absFront_0_0`
Front absortivity of construction 0 and layer 0 of a fenestration systems. As many as constructions and layers.

:file:`absBack_0_0`
Back absortivity of construction 0 and layer 0 of a fenestration systems. As many as constructions and layers.

:file:`cntrlMtx`
Matrix control file. :file:`-mtxCntrl`

:file:`shadingSch`
Schedule of control states. :file:`-schCntrl`

:file:`outside`
Radiance definition of an outdoor scene. :file:`-outside`

:file:`VARIABLES`: This section of the configuration file contains all the unidimensional inputs of the program. These are the following:

:file:`lat`
Site latitude is degrees north (use negative angle for south latitude).

:file:`lon`
Site longitude is degrees west (use negative angle for east longitude).        

:file:`tzone` (e.g. :file:`-15.00`) 
The site standard meridian is degrees west of Greenwich (use negative angle for east).      

:file:`iniMonth`
Begin Month           

:file:`iniDay`
Begin Day of Month           

:file:`endMonth`
End Month         

:file:`endDay` 
End Day of Month             

:file:`iniDayWeek`
Day of Week for first of January {1-Monday, 2-Tuesday}          

:file:`volume`
Room volume {m3}

:file:`floorArea`
Floor area {m2}

:file:`floor` (e.g. :file:`0`) 
Floor surface {id of the surface considered to be the floor, 0-first}

:file:`grndAlb`
Ground albedo

:file:`numConWin`
Number of window construction

:file:`powerLight`
Lights Watts per Zone Floor Area [W/m2]. :file:`-therm`            

:file:`airExch`
Infiltration/ventilation {ACH}. :file:`-therm`       

:file:`iniTemp`
Initial temperature (K). :file:`-therm`            

:file:`powerEquip`
Equipment Watts per Zone Floor Area {W/m2}. :file:`-therm`             

:file:`radFracEquip`
Equipment Fraction Radiant. :file:`-therm`        

:file:`powerPeople`
People Watts per Zone Floor Area {W/m2}. :file:`-therm`               

:file:`radFracPeople`
People Fraction Radiant. :file:`-therm`    

:file:`radFracLight`
Lights Fraction Radiant. :file:`-therm`

:file:`numPhotosensX` 
Grid of photosensors: number of sensors along the length of the floor surface. :file:`-grid`

:file:`numPhotosensY`
grid of photosensors: number of sensors along the width of the floor surface. :file:`-grid`

:file:`photosensHeight` 
grid of photosensors: height of sensors from Radiance coordinate reference (i.e. thickness of floor + distance from floor). :file:`-grid`

:file:`gridXOffset`
X-offset for photosensor grid (distance from the inner surface of east and west walls). :file:`-grid`

:file:`gridYOffset` 
Y-offset for photosensor grid (distance from the inner surface of south and north walls). :file:`-grid`

:file:`rotAng`
Building orientation (relative to true north {deg}, clockwise is negative)

:file:`length`
Zone length (East-West Axis [m]) (outer perimeter). :file:`-shoeBox`              

:file:`width`
Zone width (North-South Axis [m]) (outer perimeter). :file:`-shoeBox`                

:file:`height`
Ceiling Height [m] (outer perimeter). :file:`-shoeBox`           
      
:file:`albWall`
Wall albedo. :file:`-shoeBox`

:file:`albCeiling`
Ceiling albedo. :file:`-shoeBox`

:file:`albFloor`
Floor albedo. :file:`-shoeBox`

:file:`bcSouth`
Boundary condition southWall {0-interior,1-exterior}. :file:`-shoeBox`

:file:`bcEast`
Boundary condition eastWall {0-interior,1-exterior}. :file:`-shoeBox`

:file:`bcNorth`
Boundary condition northWall {0-interior,1-exterior}. :file:`-shoeBox`

:file:`bcWest`
Boundary condition westWall {0-interior,1-exterior}. :file:`-shoeBox`

:file:`bcCeiling`
Boundary condition ceiling {0-interior,1-exterior}. :file:`-shoeBox`

:file:`bcFloor`
Boundary condition floor {0-interior,1-exterior}. :file:`-shoeBox`

:file:`conSouth`
Construction ID southWall. \seealso{construction}. :file:`-shoeBox`

:file:`conEast`
Construction ID eastWall. :file:`-shoeBox`

:file:`conNorth`
Construction ID northWall. :file:`-shoeBox`

:file:`conWest`
Construction ID westWall. :file:`-shoeBox`

:file:`conCeiling`
Construction ID ceiling. :file:`-shoeBox`

:file:`conFloor`
Construction ID floor. :file:`-shoeBox`

:file:`numPanes`
Number of window panes. :file:`-klems`

:file:`numBsdfLay`
Number of BSDF layers. :file:`-klems`

:file:`conOpen`
Construction of open position. :file:`-mtxCntrl`

.. _other-input-files:

Other Input Files
=============================

.. _meteorological-data:

Meteorological data
-----------------

.. _geometry:

Geometry
-----------------

.. _material_database:

Material database
-----------------

.. _constructions:

Constructions
-----------------

.. _artificial_lighting_control:

Artificial lighting control
-----------------

.. _schedules:

Schedules
-----------------

.. _shading-control:

Shading Control
-----------------

.. _output-files:

Output files
=============================

.. _tutorial:

Steps to define a case study
=============================
