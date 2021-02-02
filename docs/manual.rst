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
Simplified DGP and Enhanced Simplified DGP methods. Thermal calculations are carried out 
according to the Kuhn2011 model for the heat transfer through the 
window, which additionally requires a U-value and Directional Solar Heat
Coefficients (DSHGC) as the thermal characterization of a fenestration 
system. Transfer functions and the energy balance method make it 
possible to calculate indoor air temperature and energy loads. Some 
auxiliary programs are also available to calculate BSDF and DSHGC for 
simplified systems composed of a shading device and a glazing unit. 

.. _examples:

Examples
=============================

For a completely new case study, where also the weather data file is 
used for the first time, a full Fener command would look like::

  > python3 ./src/fener.py {directory_path}/config.fnr -c 20 -meteo -geo -shoeBox -daylight -grid -lightSch -therm -glareSimpl

This case study consists of 20-core simulation for a shoe-box room 
geometry. The program generates the meteorological files, the Radiance 
geometry (including a grid of horizontal photosensors) and the 
three-phase method matrices. The program runs daylight, glare and 
thermal simulations. A schedule of artificial light operation is 
generated and used in the thermal simulation. 

Once the meteorological and geometry files are created, running 
variations of the same case study that share these files can be realized 
by running::

  > python3 ./src/fener.py {directory_path}/config.fnr -daylight -grid -lightSch -therm -glareSimpl

On the other hand, if the goal is to perform just a daylight simulation 
for a grid of horizontal photosensors (no thermal or glare 
calculations), the command line is the one below. In this case, 
meteorological files, room geometry and three-phase method matrices have 
already been created from previous simulations::

  > python3 ./src/fener.py {directory_path}/config.fnr -daylight -grid

.. _program_options:

Program Options
=============================

In the program commands above, :file:`fener.py` is the main routine of the program, which contains calls to all other routines, and :file:`config.fnr` is the input data file. The rest are options of the program described below:

Option :file:`-c` sets the number of cores (processes) to accelerate computation on a shared memory machine. This option is highly recommended to be used when using Radiance (e.g. generation of three-phase method matrices (option :file:`-geo`) or enhanced-DGP glare calculations (option :file:`-glare`).

Option :file:`-meteo` generates required files from the meteorological file. This option must be used when a new meteorological data file is used. Generated files are stored in the same folder as the meteorological file and are available for other simulation with the same meteorological data.

Option :file:`-geo` generates the Radiance geometry and three-phase method matrices of the case study from the program inputs. This option must be used when a new geometry (or a change in the geometry) is specified. The files generated are stored in the folder :file:`workDir/` and are available for other simulations with the same geometry.

Option :file:`-shoeBox` builds a shoe-box geometry (rotatable rectangular space). This option creates a file :file:`surfaces.dat`, which has the definition of opaque surfaces.

Option :file:`-grid` builds a horizontal grid of illuminance sensors. This option must be used combined with the options :file:`-geo` and :file:`-daylight` or with :file:`-geo` and :file:`-df`. 

Option :file:`-display` visualizes the geometry using :file:`objline`, radiance quick preview program. Simulation is paused during displaying. This option must be used combined with the option :file:`-geo`.

Option :file:`-daylight` runs a daylighting simulation.
Without the :file:`-grid` option, a Radiance description of horizontal illuminance sensors must be provided and referred in the :file:`config.fnr` file.

Option :file:`-lightSch` calculates a schedule of artificial lighting operation from the daylight simulation.

Option :file:`-glareSimpl` runs a glare simulation based on the Simplified-DGP method. The Simplified DGP is proportional to the vertical illuminance at predefined sensor positions. 
A Radiance description of vertical illuminance sensors must be provided and referred in the :file:`config.fnr` file.

Option :file:`-glare` runs a glare simulation based on the Enhanced Simplified DGP method. The glare model requires a Radiance simulation every timestep, which implies an added computational cost to the simulations. A Radiance geometry and material definition of the fenestration system must be provided and referred in the :file:`config.fnr` file.

Option :file:`-therm` runs a thermal simulation. The transient heat transfer through opaque elements (walls, ceiling) is solved based on transfer functions. The heat transfer through fenestration systems is solved based on the Black-Box model (Kuhn et al 2011).
Angular-dependent g-values or Directional Solar Heat Coefficients (DSHGC) of the fenestration system must be provided and referred in the :file:`config.fnr` file.

Option :file:`-thermComf` runs a thermal comfort calculation, based on the computer model of ASHRAE 55-2010. This module calculates a seven-point thermal sensation scale (PMV) and a percentage of dissatisfied occupants (PPD) (ANSI/ASHRAE Standard 55-2010).
The option :file:`-thermComf` can only be run together with :file:`-therm`.

Option :file:`-df` runs a daylighting factor simulation. The sky distribution for this calculation corresponds to a standard CIE overcast day. 

Option :file:`-da` performs a daylight autonomy calculation. This option must be used combined with the option :file:`-daylight`.

Option :file:`-klems` generates a BSDF dataset and angular-dependent layer absortivities of a fenestration system from the BSDF of the individual layers based on the Klems' method. BSDF files are saved according to the paths indicated in the :file:`config.fnr` file. 

Option :file:`-outside` includes an outdoor scene in Radiance format for the Daylight matrix calculation. A file containing the outdoor scene must be specified in the :file:`config.fnr` file. 

Option :file:`-schCntrl` controls a fenestration system based on a predefined schedule of window construction. The schedule must be specified in the :file:`config.fnr` file.

Option :file:`-mtxCntrl` controls a fenestration system based on a matrix of conditions that depend on internal variables (e.g. indoor temperature, solar radiation on the facade, etc.).

.. _input-data-file:

Input Data File
=============================

The :file:`config.fnr` file is the input data file of Fener. It is composed of two sections: :file:`PATHS` and :file:`VARIABLES`. It may include comments, prefixed by the character #. Unless otherwise indicated, units are always IS (e.g. m, K, W, etc.). The inputs required by each program option are indicated below.

:file:`PATHS`

This section of the configuration file tells the program  the input, output and working directories, as well as the files where the multidimensional inputs are defined. In Fener, multidimensional inputs are always defined in text files. For example, the construction of a wall is a multidimensional input variable that contains the material indexes for each layer of the wall. Multidimensional inputs are the following:

:file:`meteo` 
Meteorological information: air temperature, direct normal solar irradiation, diffuse horizontal solar irradiation, incoming longwave radiation and wind speed.

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
Calorimetric dataset of the fenestration system. 

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
Radiance geometry and material of the fenestration system. Option :file:`-glare`

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
Front absortivity of construction 0 and layer 0 of a fenestration systems. As many as constructions and layers. Option :file:`-therm`

:file:`absBack_0_0`
Back absortivity of construction 0 and layer 0 of a fenestration systems. As many as constructions and layers. Option :file:`-therm`

:file:`cntrlMtx`
Matrix control file. Option :file:`-mtxCntrl`

:file:`shadingSch`
Schedule of control states.  Option:file:`-schCntrl`

:file:`outside`
Radiance definition of an outdoor scene. Option :file:`-outside`

:file:`VARIABLES`: This section of the configuration file contains all the unidimensional inputs of the program. These are the following:

:file:`lat`
Site latitude is degrees north (use negative angle for south latitude).

:file:`lon`
Site longitude is degrees west (use negative angle for east longitude).        

:file:`tzone` (e.g. :file:`-15.00`) 
The site standard meridian is degrees west of Greenwich (use negative angle for east).      

:file:`iniMonth`
Begin Month of the simulation.

:file:`iniDay`
Begin Day of Month of the simulation.        

:file:`endMonth`
End Month of the simulation.         

:file:`endDay` 
End Day of Month of the simulation.            

:file:`iniDayWeek`
Day of Week for first of January {1-Monday, 2-Tuesday, ...}          

:file:`volume`
Room volume [m3]

:file:`floorArea`
Floor area [m2]

:file:`floor` (e.g. :file:`0`) 
Floor surface {id of the surface considered to be the floor, 0-first}

:file:`grndAlb`
Ground albedo.

:file:`numConWin`
Number of window constructions.

:file:`powerLight`
Lights Watts per Zone Floor Area [W/m2]. :file:`-therm`            

:file:`airExch`
Infiltration/ventilation [ACH]. :file:`-therm`       

:file:`iniTemp`
Initial temperature [K]. :file:`-therm`            

:file:`powerEquip`
Equipment Watts per Zone Floor Area [W/m2]. :file:`-therm`             

:file:`radFracEquip`
Equipment Fraction Radiant. :file:`-therm`        

:file:`powerPeople`
People Watts per Zone Floor Area [W/m2]. :file:`-therm`               

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
Building orientation (relative to true north [deg], clockwise is negative)

:file:`length`
Zone length (East-West Axis [m]) (indoor perimeter). :file:`-shoeBox`              

:file:`width`
Zone width (North-South Axis [m]) (indoor perimeter). :file:`-shoeBox`                

:file:`height`
Ceiling Height [m] (indoor perimeter). :file:`-shoeBox`           
      
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
Construction ID southWall. :file:`-shoeBox`

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

The meteorological data file used in Fener simulations must be in :file:`.epw` format. The format dictionary is decribed in :file:`http://apps1.eere.energy.gov/buildings/energyplus/pdfs/auxiliaryprograms.pdf`. Available meteorological files for different sites around the world can be found here:
:file:`http://apps1.eere.energy.gov/buildings/energyplus/weatherdata_about.cfm`

Whenever a new meteorological file is used, new files must be generated for the simulations through the option :file:`-meteo`. The files generated are the following: :file:`.wea` and :file:`.smx`, by using the Perez all-weather model for the visible and the solar range. A file of hourly values for the solar altitude and azimuth are also generated.

.. _geometry:

Geometry
-----------------

In Fener, the geometry information is contained in the following three files whose paths are defined in the config file:

Opaque surfaces (e.g. walls, ceiling, floor, etc.), hereafter referred as :file:`surfaces`.

Window frames, hereafter referred as :file:`frames`.

Window translucent areas, hereafter referred as :file:`windows`.

In order to understand the geometry definition of Fener, the following rules must be observed:

Each geometry input file has one header line not read by the program.
Every subsequent line of the file refers to a new element, i.e. if the window file has three lines (apart from the header), that means three windows are defined.

The surface file (:file:`surf`) is composed of the following fields: :file:`length, height, thickness, tx, ty, tz, rx (deg), ry (deg), rz(deg), ExtBoundaryCond, svf, exterior albedo, interior albedo and construction`. Each surface is built on coordinates of the origin(0,0,0) and then moved according to its translation (tx, ty, tz) and rotation (rx, ry, rz) parameters.

The frame file (:file:`frame`) is composed of the following fields: :file:`length (m), height (m), thickness (m), surface, x-offset (m), z-offset (m), out reveal (m), U-value, exterior albedo, interior albedo and svf`. The field :file:`surface` indicates the ID number of the containing surface. The fields :file:`x-offset` and :file:`z-offset` define the position of the frame with respect to the lower-left corner of the surface (from outside). The field :file:`outside reveal` defines the position of the frame with respect to the outer plane of the surface. 

The window file (:file:`win`) is composed of the following fields: :file:`length (m), height (m), frame, x-offset (m), z-offset (m), out reveal(m), svf and construction`. The field :file:`frame` indicates the ID number of the containing frame. The fields :file:`x-offset` and :file:`z-offset` define the position of the window with respect to the lower-left corner of the frame (from outside). The field :file:`outside reveal` defines the position of the window with respect to the outer plane of the surface. A window element is considered infinitely thin.  

Windows are contained in frames, and frames are contained in surfaces. Therefore, translation and rotation parameters defined for one surface also affect the frames and windows contained in that surface. 

This geometry information is used by the program to generate a Radiance geometry in the :file:`workDir` folder especified in the config file. The program also creates a network of irradiance sensor points around the surfaces (if :file:`-therm`) and a grid of horizontal illuminance sensor points (if :file:`-grid`). New three-phase method matrices are generated in the :file:`workDir` folder. The option :file:`-geo` deletes all the previous files in the :file:`workDir` folder. 

The :file:`-shoeBox` option can be used to create the file :file:`surf` as a rotatable rectangular shoe-box space. By using this options, indoor dimensions must be provided and the thickness of the enclosure is assumed to be 0.15 m. Note that the :file:`frame` and :file:`win` files must still be manually created as described above.

.. _material_database:

Material database
-----------------

Material information is user-specified in the :file:`config.fnr` file. Each material file has one header line not read by the program. The following files can be found in the database:

- Opaque material (:file:`matOpaque`).This file contains the following fields: thickness, thermal conductivity, density and specific heat (surface properties are defined in the file :file:`surf`).
- Glazing material (:file:`matGlz`). This file contains the following fields: thickness, solar transmittance, solar front reflectance (outside), solar back reflectance, visible transmittance, visible front reflectance, visible back reflectance, IR transmittance, front emissivity, back emissivity, conductivity, q-parameter (Roos model) and volumetric heat capacity. Option :file:`-klems`.
- Gas material (:file:`matGas`). This file contains the following fields: thickness, fraction of air and fraction of argon. Option :file:`-klems`.
- This file contains the following fields: thickness, thermal conductivity, front emissivity, back emissivity, top opening multiplier, bottom opening multiplier, left side opening multiplier, right side opening multiplier, front opening multiplier and volumetric heat capacity. Option :file:`-klems`.


.. _constructions:

Constructions
-----------------

Construction information is contained in the following files whose paths are defined in the config file:

Opaque surfaces (:file:`constOpaq`).
Glazing areas (:file:`constWin`).


In order to understand the construction definition of Fener, the following rules must be observed:

- Each construction file has one header line not read by the program.
- Every subsequent line of the file refers to a new construction, i.e. if the window construction file has three lines (apart from the header), it means that three window constructions are defined.
- Each construction information is constructed as a 2D-array during simulations, exclusive header line.
- The first column of a construction line is the number of layers (excluding window gaps in window constructions).
- In opaque constructions (except first column), sequence of column indicates the ID of the materials from outside to inside (from outside to inside). The ID of the materials are refered from  :file:`matOpaque.dat`.
- Constructions applied to surfaces with adiabatic boundary conditions must be symmetric with respect to the vertical center line in the cross section of the surface.
- In window constructions (except first column), the numbers after the first one come in pairs. 
The first number of the pair is the material type:
:file:`0-gas, 1-glazing, 2-shade and 4-BSDF`.
The second number of the pair is the material ID. 
For example, if the first number of the pair is 1 (=glazing), the second number of the pair indicates the ID of the material glazing (line number in :file:`matGlazing.dat`). 


.. _artificial_lighting_control:

Artificial lighting control
-----------------

Artificial lighting control information is contained in the :file:`lightCntrl` file whose path are defined in the :file:`config.fnr` file. This file has one header line and as many other lines as control photosensors. For each control photosensor, the following information must be defined:

- Control photosensor ID according to the list of photosensors defined in the :file:`config` file. For the option :file:`-grid`, [0-photosensor closest to south-west corner,1-next to east,etc.]
- Fraction of the floor area controlled by the photosensor [0-1].
- Illuminance setpoint [lux]
- Control type [0-continuous dimming control,1-stepped control]
- Control parameter [Dimming: minimum light output fraction; Stepped: number of steps].


.. _schedules:

Schedules
-----------------

Schedule information is contained in the :file:`config` file. Each schedule file contains a single column of hourly values for each hour of the year (total: 8760 values, no header). This values are multiplied by the corresponding power (if applies) defined in the :file:`config.fnr` file.

Compact schedules:

The program :file:`genSch.py` creates a schedule file in the right format from a schedule file in compact form. A schedule in compact form is composed of the following information:

- One-line header not read by the program.
- One line with the number of hour steps for working days and for weekend days (comma separated).
- As many lines as number of hour steps for working days. For each line there are two values: the end hour of the step (the start hour is the previous end hour or 00:00 for the first hour), and the schedule value.
- As many lines as number of hour steps for weekend days with the same pair of values described above. 

For example, a compact definition for heating setpoints such as::

     One-line Header
     3,1
     6,15.0
     18,20.0
     24,15.0
     24,15.0

means:

For: Working days (Mo-Fr), 3 hour steps:
(from 01:00) Until: 06:00, Heating setpoing: 15 C,
(from 07:00) Until: 18:00, Heating setpoing: 20 C,
(from 19:00) Until: 24:00, Heating setpoing: 15 C.
For weekend days (Sa-Su), 1 hour step:
(from 01:00) Until: 24:00, Heating setpoing: 15 C.

The program :file:`genSch.py` requires a schedule file in compact form and the name of the newly generated schedule file::

    > genSch.py schSpHeatConst heatSetpointConst.dat

.. _shading-control:

Shading Control
-----------------

The most simple option :file:`-schCntrl`is to assign a schedule of control states. This must be referred in the :file:`config.fnr` file and follow the format of other schedules.

The option :file:`-mtxCntrl` defines control algorithms that depend on simulation variables, such occupation, solar radiation on the facade or daylighting levels, must be defined through the file :file:`cntrlMtx` indicated in the :file:`config.fnr` file. 

The format of the :file:`cntrlMtx` file if the following:

One-line header not read by the program.
List of all the simulation variables included in the conditions. Available simulation variables are the following:
- occupation
- sunAltitude [degree]
- average workplane illuminance [lux]
- dgp
- indoor air temperature [K]
- average 24-hour temperature [K]
- exterior irradiance on south facade [W/m2]

List of setpoints for the different variables above which a condition is fulfilled. For example, for a system that has only two states (0: OFF, 1: ON) and an algorithm that is activated when there is occupation and the solar radiation on the facade is higher than 150 W/m2, we will introduce 1 under 'occupation' and 150 under 'Exterior irradiance on building facade'.
Matrix that relates conditions on simulation variables with control states. The matrix has as many rows as conditions and as many columns as variables included in the conditions. In the example before, the control matrix will have one row and two columns plus an additional column indicating the control state index:
1 1 1
For all the other situations, the control will take the default system defined with the variable :file:`conDef` in the :file:`config` file, in this case :file:`conDef = 0` (OFF). The corresponding pseudo-code is the following::

   if occupation and solar_radation > 150W/m2:
      ON
   else:
      OFF

An elaborated example of shading control algorithm is presented here. Again, we have a system with only two states (0: OPEN, 1: CLOSE). When the room is unoccupied, the algorithm compares the simulated indoor temperature with a low temperature setpoint in order to decide whether to activate the shades during the day blocking solar heat gains or to deactivate them during the night enhancing heat transfer through the window. When the room is occupied, a minimum daylighting level is imposed before closing the shade. Once the daylight condition is fulfilled, the algorithm checks the maximum vertical illuminance and the indoor air temperature. If any of these variables reaches a certain threshold, the shades are activated. The resulting algorithm is here written in pseudo-code:

if occupation:
   if average_workplane_illuminance > 400 lux:
      if indoor_air_temperature > 25C or
         max_vertical_illuminance > 3500 lux:
         CLOSE
      else: 
         OPEN
   else:
      OPEN
else:
   if night:
      if indoor_air_temperature > 19C:
         OPEN
      else:
         CLOSE
   else: 
      if indoor_air_temperature > 19C:
         CLOSE
      else:
         OPEN

The corresponding :file:`cntrlMtx` file if the following::

   # one-line header
   0,2,3,4,1,4
   1,400.,0.4,298.15,1,292.15
   0,0,0,0,0,0,1
   0,0,1,0,0,0,1
   0,0,0,1,0,0,1
   0,0,1,1,0,0,1
   0,1,0,0,0,0,1
   0,1,1,0,0,0,1
   0,1,0,1,0,0,1
   0,1,1,1,0,0,1
   0,0,0,0,1,1,1
   0,0,1,0,1,1,1
   0,0,0,1,1,1,1
   0,0,1,1,1,1,1
   0,1,0,0,1,1,1
   0,1,1,0,1,1,1
   0,1,0,1,1,1,1
   0,1,1,1,1,1,1
   1,1,1,0,0,0,1
   1,1,1,0,0,1,1
   1,1,0,1,0,0,1
   1,1,0,1,0,1,1
   1,1,1,1,0,0,1
   1,1,1,1,0,1,1
   1,1,1,0,1,0,1
   1,1,1,0,1,1,1
   1,1,0,1,1,0,1
   1,1,0,1,1,1,1
   1,1,1,1,1,0,1
   1,1,1,1,1,1,1

.. _output-files:

Output files
=============================

Output files are saved in the folder :file:`output`. One output file is generated per variable for the simulation period.  Output files do not contain header or time indication. Current output files are the following: 

:file:`ill.out` - [lux] Illuminance map at grid of photosensors for each simulation hours.

If option :file:`-glare`: 

:file:`illVert.out` - [lux] Vertical illuminance at glare sensor points. 

:file:`dgp.out` - [0-1]  Daylight glare probability (DGP) index at glare sensor points. 

If option :file:`-therm`:

:file:`solTrans.out` - [W] Solar irradiation transmitted through windows.

:file:`effGValue.out` - [-] Effective g-value.

:file:`convIntHeatFlx.out` - [W] Convective heat flux from internal heat gains.

:file:`radIntHeatFlx.out` - [W] Radiant heat flux from internal heat gains.

:file:`irrSurfExt.out` - [W m-2] Solar irradiation absorbed by outdoor opaque surfaces.

:file:`irrSurfInt.out` - [W m-2] Solar irradiation absorbed by indoor opaque surfaces.

:file:`infrSurfInt.out` - [W Infrared heat flux on indoor opaque surfaces.

:file:`convSurfInt.out` - [W] Convective heat flux on indoor opaque surfaces.

:file:`heatSurfInt.out` - [W] Heat flux from internal heat fluxes on indoor opaque surfaces.

:file:`tempSurfInt.out` - [K] Temperature of indoor opaque surfaces.

:file:`tempSurfExt.out` - [K] Temperature of outdoor opaque surfaces.

:file:`irrFrameExt.out` - [W m-2] Solar irradiation on outdoor frame surfaces.

:file:`irrFrameInt.out` - [W m-2] Solar irradiation on indoor frame surfaces.

:file:`irrWinAbsExt.out` - [W m-2] Solar irradiation absorbed by the outdoor layer of windows.

:file:`irrWinAbsInt.out` - [W m-2] Solar irradiation absorbed by the indoor layer of windows.

:file:`infrWinInt.out` - [W] Infrared heat flux on indoor glazing surfaces.

:file:`convWinInt.out` - [W] Convective heat flux on indoor glazing surfaces.

:file:`heatWinInt.out` - [W] Heat flux from internal heat fluxes on indoor glazing surfaces.

:file:`tempWinInt.out` - [K] Temperature of indoor glazing surfaces.

:file:`tempWinExt.out` - [K] Temperature of outdoor glazing surfaces.

:file:`tempInAir.out` - [K] Indoor air temperature.

:file:`energyDemand.out` - [W] Zone energy demand (cooling negative).

If option :file:`-thermComf`:

:file:`pmv.out` - [-] seven-point thermal sensation scale (PMV) of five equally distributed sensors in the room (ASHRAE 55-2010).

:file:`ppd.out` - [%] percentage of thermal dissatisfied occupants (PPD) of five equally distributed sensors in the room (ASHRAE 55-2010).

:file:`tWalls.out` - [K] Temperature of combined opaque and glazing surfaces. 

:file:`tRad.out` - [K] Mean Radiative Temperature of five equally distributed sensors in the room. 

:file:`RHOut.out` - [%] Relative Humidity of the room exterior.

:file:`RHIn.out` - [%] Relative Humidity of the room interior.

If option :file:`-da`:

:file:`da.out` - [%] Spatial daylight autonomy at photosensors.

If option :file:`-df`:

:file:`df.out` - [%] Daylight factor at photosensors.

:file:`-schCntrl` :file:`-mtxCntrl` :

:file:`conIndex.out` - Construction schedule of glazing areas.

If option :file:`-lightSch`:

:file:`lightSch.out` - Schedule of artificial lighting for simulation period.
 
:file:`lights.dat` - Schedule of artificial lighting for the total hours of the year to be used by EnergyPlus. This file is the same as the input file :file:`lightSch` in case the option :file:`-lightSch` is not used.

.. _tutorial:

Steps to define a case study
=============================

The steps to define a new case study and run a coupled daylighting and thermal simulation are the following:

- Copy the required meteorological data file in your input folder.
- Modify the input data in :file:`config.fnr` file.
- Modify the artificial lighting control file :file:`lightControl.dat`.
- Modify the geometry files :file:`frame` and :file:`win`.
- Modify the construction file :file:`constOpaq`.
- Check if the required materials are defined in the material files.
- Modify the schedule compact files.
- Run the program :file:`genSch.py`::

    > genSch.py schSpHeat heatSetpoint.dat
    
- Run the program::

    > fener.py config -c 20 -meteo -shoeBox -grid -geo -daylight -lightSch -therm
