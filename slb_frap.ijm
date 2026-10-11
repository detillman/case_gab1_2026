// for each file, auto-adjust contrast for first frap frame
// draw one circular roi around bleached area
// rename each roi with file basename
// save all rois in a single RoiSet.zip file

// defining roi_finder function
function roi_finder(file_name) { 
	number_rois = roiManager("Count"); 
	for (roi_index=0; roi_index<number_rois; roi_index++) { 
		roiManager("Select", roi_index); 
		roi_name = Roi.getName(); 
		if (matches(file_name, roi_name)) { 
			return roi_index; 
		} 
	} 
	return -1; 
} 

// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// set parameters
run("Set Measurements...", "area mean min redirect=None decimal=4");

// open and sort rois
open(dir + "RoiSet.zip")
roiManager("Sort");

// run macro in batch mode
setBatchMode(true);
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	strip_name = replace(name, ".tif", "");
	if (endsWith(path, "tif")) {

		// open and auto-ajdust stack
		open(strip_name + ".tif");
		rename("in_use");
		setSlice(1);
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		
		// get stack properties		
		frame_interval = Stack.getFrameInterval();
		number_slices = nSlices();
		
		// make empty arrays
		name_results = newArray(number_slices);
		slice_results = newArray(number_slices);
		time_results = newArray(number_slices);
		totalMean_results = newArray(number_slices);
		totalArea_results = newArray(number_slices);
		roiMean_results = newArray(number_slices);
		roiArea_results = newArray(number_slices);
		
		// measure total signal across all slices
		for (slice=1; slice<(number_slices+1); slice++) {
			setSlice(slice);
			run("Measure");
			
			// store total values in arrays
			name_results[slice-1] = strip_name;
			slice_results[slice-1] = slice;
			time_results[slice-1] = (slice-1)*frame_interval;
			totalMean_results[slice-1] = getResult("Mean", slice-1);
			totalArea_results[slice-1] = getResult("Area", slice-1);
		}
		close("Results");
		
		// find roi that matches file name
		roi_index = roi_finder(strip_name);
		roiManager("Select", roi_index);
		
		// measure roi signal across all slices
		for (slice=1; slice<(number_slices+1); slice++) {
			setSlice(slice);
			run("Measure");
			
			// store roi values in arrays
			roiMean_results[slice-1] = getResult("Mean", slice-1);
			roiArea_results[slice-1] = getResult("Area", slice-1);
		}
		close("Results");
		
		// store arrays in table and save
		Table.create("saveTable");
		Table.setColumn("file", name_results);
		Table.setColumn("slice", slice_results);
		Table.setColumn("time", time_results);
		Table.setColumn("totalMean", totalMean_results);
		Table.setColumn("totalArea", totalArea_results);
		Table.setColumn("roiMean", roiMean_results);
		Table.setColumn("roiArea", roiArea_results);
		Table.save(dir + "/" + strip_name + ".csv");
		close("saveTable");
		close("in_use");
	}
}

// close ROI manager and end
close("ROI Manager");
beep();