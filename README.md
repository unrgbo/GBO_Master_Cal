Series of routines for creating master calibration frames and searching calibration files.
	
Many of the programs define paths in a format for windows; however, they often have alternate linux/unix style paths commented out. Just update all paths with the desired format and names for your directory tree. There is a lot of commenting that still needs to be added.

master_cal_weather_v2.py: is the primary script for searching a defined directory and finding all available calibration frames, then produces master frames for sets. It is currently setup very generalized and produces masters that are not bias or dark corrected. There is options to do this and create true bias and dark corrected normalized flats; however, in order to keep them usable in any software, they currently produce masters that are simply median combined images that still contain bias and dark current. This can easily be adjusted within each routine [master_bias(), master_dark(), master_flat()].

calgui.py: is a useful script for searching calibration master frames based on master_cal_weather_v2.py output for calibration file info. It is currently set up as a GUI but can be ran from the terminal with a couple of quick changes.
