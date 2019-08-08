# FMD-disease
Project name: Benefits of considering landscape context in disease control

Description: The project simulates the evolution of a World of animals and the spread of a disease in this World at each time step (one day in our simulation). Animals
are distributed in farms. We propose different strategies to tackle the spread of the disease and to minimize the number of killed animals.

Contents : 
The pdf file is the article, the .py file contains Python file. The .csv file is a data set taken from Sellman's article, we import it 
in the Python file.
The .gif shows an animation of the initial world during 10 days.

Installation :
We use the disease parameters from the Sellman's article. 
The project is semi-oriented object : the class farm has 4 parameters and some functions.
To generate the World, we can use the function Random_World (generate farms with randon positions) or use the main text to construct a World from a Data file .csv.
Iteration_0 updates the world with any management. Iteration_1, Iteration_2, Iteration_3 are iteration functions of the different strategies.
Draw_Susceptibles, Draw_Infectious, Draw_Culled and Draw_World can be used to draw the world.
Iteration_2 uses the functions CloserSusceptibles, Simulate_World and Rank_Susceptibles.
Several_Worlds and ListToArray are auxiliary functions that help for the plots of the different strategies. 

Credits :
Sellman's article was the previous work for FMD disease : https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006086

License :
MIT license
