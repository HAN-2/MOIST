# MOIST--a Matlab based One-dimensional Isotope and Soil water Transport model
(The update of this guide will be slow because I have limited time to write in detail. However, I will try my best to finish this document.)


Codes for theoretical tests, semi-analytical tests, and lysimeter validations. They are only guarented to regenerate the results in the manuscript ('Title of the manuscript will appear here').<br>
If you have any further questions, welcome to send an email to Han (haf033@usask.ca Please use MOIST in the subject and thank you for your interest!).

Validation data can be download at: <br>
EPFL:https://zenodo.org/record/4037240#.Y029l3bMKUk<br>
HBLFA:https://www.pc-progress.com/en/Default.aspx?h1d-lib-isotope


# HOW TO RUN THESE CODES:
Each scenario has a 'Main_program' or 'Main_program_2' file. It is runable when all the function files are located under the same path.<br>
(If there are both 'Main_profram' and 'Main_program_2', 'Main_program_2' is preferred)

# MAIN OUTPUT VARIABLES:<br>
There are plenty of output variables after the calculation is done. Following variables could be considered the most important. <br>
More explainations about variables will be updated. <br>

### UPDATING.... ###

z_theta: output of soil water content profiles<br>

z_cil:   output of isotope profiles<br>

z_T:     output of temperature profiles<br>

z_h:     output of soil water head profiles<br>

z_qevap: evaporation flux matrix<br>

z_dt:    time step matrix<br>

for semi-analytical tests:<br>
zBA1: results of isotope profiles under isotherm and saturated conditions<br>

zBA:  results of isotope profiles under non-isotherm and non-saturated conditions<br>

# HOW TO PREPARE DATA:
All needed data are integrated into the EXCEL file, like 'Magali.xlsx' and 'stumpp.xlsx' in short and long term validation folders, respectively. <br>
There are four sections:<br>
## 'Rain_record':  <br>
A and B column are the start and end points of each time interal (in second), repectively.<br>
C column is the rainfall amount within each time interval.<br>

## 'Daily_climate_record':<br>
this section is usually arranged as:<br>
A -->Month; B -->Day; C -->Year<br>	
D -->Tmax; E -->Tmin; F -->Tavg<br>	
G -->sunshinehour(h); H -->Daily_rh; I -->u (m/s)<br>	
J -->solar  radiation(W/m2)<br>

Note that the daily climate input also depends on data availbility, it is not always managed as above. You can customize the 'climate_function' according to your own needs (start from 'if X5ip2.day_flag'). The length of each column should be equal. If there are missing data, intrepolation can be considered.

## ???initial_condition???<br>
A -->Initial soil water content (m3/m3); B -->Initial temperature profile <br>	
C --> emplty column; D --> Initial soil water isotoic compositions <br>	

The length of each column should be equal and consistent to the intervals that assigned from main program.

## 'Rainfall_isotope'<br>
A --> Rainfall amount (m/s)<br>	
B --> isotopic compositions (???) <br>	

# WHEN THE CODE IS RUNNING:
When the code is running, there is a dynamic number (%) at the left corner of command window, which is helpful to check if the program is stucked or the roughly running time. 

# HOW TO PLOT:
There are plotting commands at the end of main program of the short term validation. Also, I will post the commands here to illustrate how to extract temporal and spatial soil water and isotope data for further visualization.
