
# Data Analysis for Folding@Home

## Authors
Samir Singh, Vincent Voelz

### Copyright
Voelz Lab, Temple University


# Description

This project contains several scripts that are used for analyzing data on the Folding@Home server and for now caters specifically to `vav3`. These  
scripts are meant to be modular and have classes and objects added to them for any further analysis that is needed.

## analysis\_featurizer.py

`analysis_featurizer.py` is the base scripts for this whole project. It contains all the classes needed for analyzing data.

### Featurizer Class

The Featurizer class is the parent class of all classes in the script. You supply it a project number, run number, clone number, and gen number so it can analyze data. It will automatically retrieve all necessary files needed for analyzing the project. It will automatically remove periodic boundary conditions from the traj\_comp.xtc file and load in a trajectory with that and a compatible gro file found in the projects directory of vav3.

You will usually not be using the Featurizer class though. The children of this class are what is used for tailoring what type of data you want calculated. For this class to work, there will need to be a features folder in the projects directory for the project you want to be analyzed. This folder must contain an xtc.ndx file created by the gromacs make\_ndx command for selecting the index group. There also needs to be a text file called index which just needs two numbers on separate lines that selects the atom groups when the class runs the gromacs trjconv command. There also needs to be a gro file in the projects directory with solvent removed so it will be compatible with the xtc file.

#### DistanceCalculator Class

This is a class used to calculate distances between atoms in a trajectory. You can specify the description of the atoms you want the distances calculated for or you can create a features.npy file in the data directory of the project that contains all the atoms you want distances calculated for. These atom distances will automatically be saved in a distances.npy file in a folder called features in the data directory of the project. This file is a numpy array of shape run, clone, gen of the specified project. Each index holds all the distances calculated for that gen.

##### Example

```
from analysis_featurizer import DistanceCalculator
t = DistanceCalculator(14137, 0, 0, 0) # (Project Number, Run Number, Clone Number, Gen Number)
t.calculate_distances() # checks for features.npy. If it doesn't exist will use alpha and beta carbons
```

Output:
```
Loading in atom indices...
Calculating distances...
Distances saved for index [0][0][0]
array([[ 4.79062223,  4.44829369,  4.61328983, ...,  4.33823109,
         4.42707348,  4.37204933],
       [ 3.93415928,  4.02548361,  4.2844348 , ...,  3.44840217,
         3.24982047,  3.06404471],
       [ 3.8231647 ,  3.92661381,  4.22720194, ...,  3.56349421,
         3.45095229,  3.13069725],
       ..., 
       [ 6.38097954,  5.9705739 ,  6.02991724, ...,  6.17971945,
         6.12991858,  5.59185696],
       [ 6.24494314,  6.21516132,  6.30268717, ...,  6.02589273,
         6.10903883,  5.98597479],
       [ 6.1682024 ,  5.87098837,  5.8595109 , ...,  5.53846836,
         5.66585827,  5.26477337]], dtype=float32)
```

#### DihedralCalculator Class

This class is used to calculate the dihedral angles of a trajectory. In addition to specifying the specific gen inside the project, you also specify  
which angle you want calculated (omega, phi, or psi). These angles are also saved to a numpy file in the features folder of the project.

#### Example

```
from analysis_featurizer import DihedralCalculator
t = DihedralCalculator(14137, 0, 0, 0, 'phi')
t.calculate_dihedrals()
```

Output:
```
calculating dihedrals for angle phi...
Dihedrals saved for index [0][0][0]
(array([[  10,   12,   14,   20],
        [  20,   22,   24,   42],
        [  42,   44,   46,   61],
        ..., 
        [4348, 4350, 4352, 4365],
        [4365, 4367, 4369, 4377],
        [4377, 4379, 4381, 4392]]), # Dihedral Indices
 array([[-1.18036175, -1.16861141, -1.37235546, ..., -2.29823136,
         -1.26390374, -1.75710154],
        [-0.95376933, -1.28588843, -1.49756622, ..., -2.18892026,
         -1.88741553, -1.90128243],
        [-1.30906892, -1.76275051, -1.50172615, ..., -2.06015015,
         -1.46920514, -1.69564986],
        ..., 
        [-1.22993827, -1.50562763, -1.82929516, ..., -2.33535004,
         -1.76825809, -1.59249353],
        [-1.86078608, -1.51037407, -1.82607496, ..., -2.09267998,
         -1.80721879, -2.23773623],
        [-1.47531056, -1.42999482, -1.32822394, ..., -1.37763286,
         -0.99655151, -2.49673223]], dtype=float32)) # Angles
```

#### COMCalculator Class

This class is used for calculating the distance between the center of mass of a group of atoms and a specified atom index. By default, the atom 
indices used for this class will come from a user created text file called 'US' in the projects directory called features which is fromatted as such:
```
[atom_index, range(atom_indices)]
```
These distances will be saved in a features folder in the data directory as a numpy file.

#### Example

```
from analysis_featurizer import COMCalculator
t = COMCalculator(14137, 0, 50, 5)
t.COM_atom_distances_for_gen()
```

Output:

```
Calculating COM for frame 0 of gen 5...
Calculating COM for frame 1 of gen 5...
Calculating COM for frame 2 of gen 5...
Calculating COM for frame 3 of gen 5...
...
Calculating COM for frame 47 of gen 5...
Calculating COM for frame 48 of gen 5...
Calculating COM for frame 49 of gen 5...
Calculating COM for frame 50 of gen 5...
COM saved to index [0][50][5]

[2.9207594578724856,  # Array of distances
 2.9747724469318149,
 2.81373559030986,
 3.0265772150001391,
 ...
 2.8673497864440338,
 2.8630320181035365,
 3.0261880831602279,
 2.6778575521785379]
```

#### RMSDCalculator Class

This class is used for calculating RMSD and by default calculates RMSD with the first frame as a reference, but there is an option to supply  
your own reference structure file. The data calculated from this is saved in the data directory in the features folder.

##### Example

```
from analysis_featurizer import RMSDCalculator
t = RMSDCalculator(14188, 0, 0, 0)
t.calculate_RMSD()
```

Output:
```
Loading in atom indices...
Calculating RMSD...
Saved RMSD to index [0][0][0]

array([ 0.        ,  0.10858581,  0.12854657,  0.14666103,  0.15992635,
        0.16884001,  0.17498644,  0.1948081 ,  0.18852335,  0.17774747,
        0.18971205,  0.19377165,  0.21642178,  0.22514133,  0.57685119,
        0.22979626,  0.23111489,  0.23691319,  0.22335187,  0.22371876,
        0.22812681,  0.23151061,  0.21934234,  0.23011787,  0.22346987,
        0.22883829], dtype=float32)
```

### Analysis Class

This class is a child of the Featurizer class but is also a base class. Its use is to analyze data, such as plotting, calculated from the Calculator  
classes. For example, the RMSDAnalysis class will plot what data has already been calculated by the RMSDCalculator class. 

#### COMAnalysis Class

This class is a child of the Analysis class that is used to plot a line graph of all the Center of Mass distances that were calculated from  
COMCalculator for a clone of data. This graph will be saved in a plots folder in the data directory of the project.

##### Example

```
from analysis_featurizer import COMAnalysis
t = COMAnalysis(14137, 0, 0, 0)
t.plot_COM_atom_distances()
```

Output:
```
Plot made for run 0, clone 0 of project 14137
# Plot image TBA
```

#### RMSDAnalysis

This class is a child of the Analysis class that plots all RMSD that has been calculated by the RMSDCalculator class. The plot is a histogram and  
the default number of bins is 500, but the number of bins can be specified in the argument of the function if you want more or less bins.

##### Example

```
from analysis_featurizer import RMSDAnalysis
t = RMSDAnalysis(14137, 0, 0, 0)
t.plot_RMSD(100)
```

Output:
```
Plotted all RMSD values available for this project
# Plot image TBA
```

#### TICA Class

This class is a child of the Analysis that calculates tICA coordinates and components and creates a tICA plot for the whole project. Everything is  
saved to a tICA folder in the data directory for the project. This class is still a work in progress.

## data\_analysis.py

This script uses analysis\_featurizer.py and calculates data for every single run, clone, gen in the specified project. It will skip over any data  
that has already been calculated for the project.

##### Example

```
from data_analysis import *
distance_calc_proj(14137)
```

Output:
```
Loading in atom indices...
Calculating distances...
Distances saved for index [0][0][1]
Loading in atom indices...
Calculating distances...
Distances saved for index [0][0][2]
...
```

## readlog.py

This is a Daemon scripts which can be run in the background of the server using the following command:
```
sudo start readlog
```
This script reads in the fah-work.log and extracts what data has come back. It then retrieves the new data and calculates the atom distances using  
the  DistanceCalculator class from analysis\_featurizer.py. readlog.py has its own log file called readlog.log where you can see what is going on  
while it is running. Future plans may include expanding readlog.py to calculate other types of data or create new scripts that will calculate data  
for a project every once in a while.

##### readlog.log output

```
2019-08-23 11:16:32: The file distances.npy has used 14mb/4000mb (0.37%) of the max user specified file size
2019-08-23 11:16:32: Calculating features and saving data (14137, 36, 212, 5)...
2019-08-23 11:16:35: Saved data for project 14137 to index [36][212][5] of distances.npy
2019-08-23 11:16:50: [Errno 2] No such file or directory: '/home/server/server2/projects/Gromacs/p14153/features/index' There is most likely no ndx.gro and index file for this project (14153, 7, 128, 162)
2019-08-23 11:16:53: [Errno 2] No such file or directory: '/home/server/server2/projects/Gromacs/p14153/features/index' There is most likely no ndx.gro and index file for this project (14153, 17, 136, 136)
2019-08-23 11:17:12: The file distances.npy has used 14mb/4000mb (0.37%) of the max user specified file size
2019-08-23 11:17:12: Calculating features and saving data (14137, 13, 205, 3)...
2019-08-23 11:17:15: Saved data for project 14137 to index [13][205][3] of distances.npy
2019-08-23 11:17:16: [Errno 2] No such file or directory: '/home/server/server2/projects/Gromacs/p14153/features/index' There is most likely no ndx.gro and index file for this project (14153, 6, 112, 131)
```


































