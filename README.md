# SUTRASET is extended USGS SUTRA code considering seepage(S), evaporation(E) and tidal forces(T)

The specific processes that the code considers on top of SUTRA are:

  1. evaporation takes water away with salt remain in the soils
  2. the reduction of evaporation due to water desaturation on the surface, by considering surface resistance.
  3. the reduction of evaporation due to increase to solute concentration on the surface, according to Kelvin Equation.
  4. the reduction of evaporation due to salt precipitation on the surface as crust, by considering salt resistance
  5. the removal of residual liquid water due to evaporation (by extending SWCC and Kr function below residual liquid water saturation level)
  6. a dynamic soil surface boundary conditions with three phases (1. Hydrostatic pressure boundary for water submerged surfaces, 2. Seepage for saturated soil surface above water level, 3 evaporation for unsaturated soil surface above water level). These three phases are determined using non-linear iterations to avoid phase lag (water spraying truck effect!).

Vapor flow in porous media is not considered in SUTRASET

  
Examples:
the result of a 1-D calibration case considering evaporation from a soil column with a saline water table at the bottom:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/Ny3e0dCjsOM/0.jpg)](https://www.youtube.com/watch?v=Ny3e0dCjsOM)

The case is available at examples/Fujimaki/MassarLoamy and /examples/Fujimaki/Toyoura


The mathematical explanation is given in the paper:
Shen, C., Zhang, C., Xin, P., Kong, J., & Li, L. (2018). Salt Dynamics in Coastal Marshes: Formation of Hypersaline Zones. Water Resources Research, 1–18. https://doi.org/10.1029/2017WR022021


The result of a 2-D case considering tide, evaporation in coastal wetland

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/y01Bo0dyTFE/0.jpg)](https://www.youtube.com/watch?v=y01Bo0dyTFE)

the case is available at examples/evaporation_salt_marsh/sandyloam_evt4


## How to use the code:

#### 1. Under UNIX environment where gfortran, git are available.

Issue the following command

```bash
git clone https://github.com/chenmingzhang/sutraset
cd sutraset/src
make 
make all
```

You will have sutraset_gf compiled in the src folder

#### 2. If you just wish to run pre-compiled SUTRASET in WINDOWS environment:

SUTRASET has been pre-compiled as **sutraset_gf_windows.exe** in the directory *bin*. The associated **dlls** are also included in **bin** folder. This exe file has a alias (or shortcut) in folder *bin* named *sutraset.bat* that can be copied to any of your simulation folders.
(We do not suggest copying the exe (sutraset_gf_windows.exe) and associated dlls file into the simulation directories for running SUTRASET.)

To allow *sutraset.bat* to be copied and executed in any of your simulation folders, you will need to add the absolute address of directory *bin* into the system environment variables.

For example: 

If **sutraset** is located in *d:\sutraset*, the absolute address of directory *bin* is *d:\sutraset\bin*. You will need to add *d:\sutraset\bin* in to your system environment variable.

The way to add path into system environment variables can be found from here:

https://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10

After this, you will be able to copy *sutraset.bat* into any of your simulation folders and run the simulation by double clicking the file.

#### 3. If you wish to compile SUTRAET in WINDOWS, with minGW available:

 After downloading the package in zip file, create a file in folder */src* called *gversion.f90* with following content

```fortran
SUBROUTINE GVERSION (K3)
  WRITE (K3,*) "this is not git controlled "
      RETURN
END SUBROUTINE
```

Then issue the command below in the folder src (with mingw) 

```bash
gfortran *.f90 *.f -O3 -o sutraset_gf_windows.exe

```

## Input requirements:
 SUTRASET treats evaporation as a sink term. The node numbers needs to be listed DATASET 17. As the change of evaporation rate due to surface saturation, concentration and solid salt depth is implemented in subroutine BCTIME, these node number will be specified as negative values. 
 Since the original SUTRA does not automatic calculate the effective surface area associated with the flux in source and sink term (dataset 17), User needs to input a volumetric inflow/outflow (m3/s) by multiplying the in/out flux (m/s) with the injecting area (m2),  Similarly, SUTRASET requires user to input the node geometry where evaporation is implemented. This is done by to specifying the surface area and the depth of the cell where evaporation is applied, after the negative node number. Please referred to example **examples/evaproation_salt_marsh/sandyloam_evt4/evapotranspiration_sandyloam.inp**
 The soil water retention parameters, air conditions, and resistance parameters are listed in input file ET.inp. See the same example mentioned above.


## Rule of Thumb to help converge your solution:
1. If negative solute concentration is present, particularly below the cell where evaporation is implemented, try to reduce the diffusivity, or refine the mesh. As evaporation may lift the salinity up to solubility limit (~265 PPT for NaCl), generating large concentration gradient near the surface. This could be relieved by quickly diffuse the hypersaline water through diffusion, or add more cells to describe the details

2. non-linear iteration must be enabled, when seepage is considered.

3. if tide, seepage and evaporation is considered on the top surface, (Refer to example */examples/evaporation_salt_marsh* ) and simulation is crashed with P equal NaN, trying to adjust GNUP to minimise the recharge in to surface unsaturated zone. This phenomenon occurs in particular when soil has high permeability.  As evaporation may de-saturate the soil surface with a saturation level of 0.01, the corresponding pressure could be -1e7 Pa according to soil water retention curve. Once tides comes with PBC of 0 Pa, large pressure deficit ( in this case 1e7) could cause significant recharge (e.g., flow_in=GNUP\*(PBC-P)=0.01\*1e7=1e5m3/kg) and so blow the surface cell.
