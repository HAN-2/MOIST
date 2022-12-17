# MOIST--a Matlab based One-dimensional Isotope and Soil water Transport model
(Note: the update of this guide will be slow because I do not have time to write it in detail. However, I will try my best to explain the variables and finish this instructions. If you have any questions, you are welcome to send an email to haf033@usask.ca and please use MOIST in the subject. Thank you for your interest!)


Codes for theoretical tests, semi-analytical tests, and lysimeter validations<br>

They are only guarented to regenerate the results in the manuscript. If you have any further questions, welcome to send an email to Han.

Validation data<br>
EPFL:https://zenodo.org/record/4037240#.Y029l3bMKUk<br>
HBLFA:https://www.pc-progress.com/en/Default.aspx?h1d-lib-isotope


# HOW TO RUN THESE CODES:
Each scenario has a 'Main_program' or 'Main_program_2' file. It is runable when all the function files are located in the same path.<br>
(If there are both 'Main_profram' and 'Main_program_2', 'Main_program_2' is preferred)

# MAIN OUTPUT VARIABLES:<br>
z_theta: output of soil water content profiles<br>
z_cil: output of isotope profiles<br>
z_T: output of temperature profiles<br>
z_h: output of soil water head profiles<br>
z_qevap: evaporation flux matrix<br>
z_dt: time step matrix<br>



for semi-analytical tests:<br>
zBA1: results of isotope profiles under isotherm and saturated conditions<br>
zBA: results of isotope profiles under non-isotherm and non-saturated conditions<br>


