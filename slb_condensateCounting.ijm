// save individual pre-frap channel(s) from .lif as "_raw.tif"

// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// set parameters and enter batch mode
run("Set Measurements...", "area mean standard min redirect=None decimal=4");
setBatchMode(true);

// run macro
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	strip_name = replace(name, "_raw.tif", "");
	if (endsWith(path, "_raw.tif")) {
	
		// open tif
		open(strip_name + "_raw.tif");
		rename("working_image");
		resetMinAndMax;
		
		// get whole-fov measurements
		run("Measure");
		fov_area = getResult("Area", 0);
		fov_sd = getResult("StdDev", 0);
		fov_min = getResult("Min", 0);
		fov_max = getResult("Max", 0);
		
		// get lut values to determine channel
		getLut(reds, greens, blues);
		Array.getStatistics(greens, green_min, green_max, green_mean, green_stdDev);
		
		// calculate and apply threshold
		//thresh = fov_min + fov_sd * 5; // for both channels
		//if (green_max > 0) {
			//thresh = fov_min + fov_sd * 12; // for green
			//}
		//if (green_max == 0) {
			//thresh = fov_min + fov_sd * 5; // for magenta
		//}
		
		// or use manual threshold
		//thresh = fov_min + fov_sd * 5; // for slb domains
		if (green_max > 0) {
			thresh = 180; // for green
			}
		if (green_max == 0) {
			thresh = fov_min + fov_sd * 5; // for magenta
		}
		
		selectImage("working_image");
		setMinAndMax(fov_min, thresh);
		
		// save tif for making figures
		saveAs("Tiff", dir + "/" + strip_name + "_adjusted.tif");
		rename("working_image");
		
		// binarize with threshold and analyze particles
		setThreshold(thresh, fov_max, "raw");
		run("Analyze Particles...", "size=0.10-Infinity display add");
		
		// save rois if rois were identified
		roi_count = roiManager("count");
		if (roi_count > 0) {
			roiManager("Save", dir + "/" + strip_name + "_rois.zip");
		}
		roiManager("reset")
		close("ROI Manager");
		
		// save results
		saveAs("Results", dir + "/" + strip_name + "_results.csv");
		close("Results");
		close("*");
	}
}

beep();