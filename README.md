# bmmb_tracking_software
This Matlab software package contains the necessary steps to idenfity and track
nuclei and Golgi bodies and then combine their identification lists to link
nuclei and Golgi bodies that belong to the same cell. This information is also 
listed in nuclei_golgi_track.m.  

Steps: 
1. Image import from tiff stacks. A user interface will appear prompting
you to select two tiff video files of tracked nuclei and Golgi bodies (in
seperate channels). If your videos are not in this format you can convert
it to a tiff stack in ImageJ or modify this Matlab code directly to suit
your needs. 
2. Image pre-processing. Images can often contain backround noise or have
insufficient brightness and contrast for objects to be picked up by the
tracking code. Please modify your inputs for this section until you see
images where you can identify objects by eye. 
3. Nuclei tracking. We use the ACTIVE nuclei tracking software package
developed by the Henderson lab. 
4. Golgi tracking. This code identifies and tracks Golgi bodies and
constructs an identification list similar to the ACTIVE package. 
5. Combination of ID lists. We look at cell positions of both nuclei and
Golgi and pair them appropriately to create one master ID list. 

Following this analysis: 
The primary tool we used to look at mean and standard deviation of
orientations was the AngleSpread2 function with the following syntax: 

[FinalStd, MaxAngle]  = AngleSpread2(theta);

where theta is a vector list of orientations. This function can be found in
the Nuclear_Alignment_Video_Analysis folder. As a note to the user this
code is more to provide the initial imaging and tracking software necessary
to identify and track irregularly shaped objects such as the Golgi body.
Any futher analysis is largely based on the user's needs. 

Developed by Giuseppe Passucci from 2013 to 2017 as part of the Manning
Group at Syracuse University in collaboration with the Henderson (SU) and
Turner (SUNY Upstate) groups. 
