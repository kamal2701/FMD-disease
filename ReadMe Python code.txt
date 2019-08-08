The Python file contains the code used for FMD project.

We use the disease parameters from the Sellman's article. 
The project is semi-oriented object : the class farm has 4 parameters and some functions.
To generate the World, we can use the function Random_World (generate farms with randon positions) or use the main text to construct a World from a Data file .csv.
Iteration_0 updates the world with any management. Iteration_1, Iteration_2, Iteration_3 are iteration functions of the different strategies.
Draw_Susceptibles, Draw_Infectious, Draw_Culled and Draw_World can be used to draw the world.
Iteration_2 uses the functions CloserSusceptibles, Simulate_World and Rank_Susceptibles.
Several_Worlds and ListToArray are auxiliary functions that help for the plots of the different strategies. 