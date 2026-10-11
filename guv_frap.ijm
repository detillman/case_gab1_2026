//// don't denoise for frap analysis but consider it for representative images
//// open .tf2 file for each frap experiment
//// auto-adjust contrast on first frame
//// set lookup table to appropriate color
//// save as .tif
//// duplicate and crop guv of interest (rectangle that's twice as large as guv)
//// auto-adjust contrast
//// save as "_frap#.tif"

//// for dynamic roi
/// draw polygon roi around bleached area in last pre-bleach frame
/// deselect roi and use selection rotator to adjust roi and add to manager for each frame
/// save rois as "_frap#.zip"
/// repeat for all frapped rois in each image

//// for static roi
/// draw polygon roi around bleached area in last pre-bleach frame
/// save roi as "_frap#.roi"
/// repeat for all frapped rois in each image

//// run macro
// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// set parameters
run("Set Measurements...", "area mean min redirect=None decimal=4");

// run macro in non-batch mode to avoid roi manager problems even though it takes longer to run
setBatchMode(false);
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	
	// process dynamic roi
	if (endsWith(path, "zip")) {

		// open and auto-ajdust stack
		strip_name = replace(name, ".zip", "");
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
		
		// open, sort, and count rois
		open(strip_name + ".zip");
		roiManager("Sort");
		roi_count = roiManager("count");
		
		// measure signal in each frame-specific roi
		for (roi=0; roi<roi_count; roi++) {
			
			roiManager("Select", roi);
			
			// break if slice and roi mismatch
			current_slice = getSliceNumber();
			if (roi != current_slice - 1) {
				exit("double-check " + strip_name + "; roi = " + roi + ", but slice = " + current_slice);
			}
			
			// measure if slice and roi match
			run("Measure");
			
			// store roi values in arrays
			roiMean_results[roi] = getResult("Mean", roi);
			roiArea_results[roi] = getResult("Area", roi);
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
		
		// close everything
		close("saveTable");
		close("in_use");
		close("ROI Manager");
	}
	
	// process static roi
	if (endsWith(path, "roi")) {
		
		// open and auto-ajdust stack and roi
		strip_name = replace(name, ".roi", "");
		open(strip_name + ".tif");
		rename("in_use");
		setSlice(1);
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		
		// open roi and add to manager
		open(strip_name + ".roi");
		roiManager("Add");
		
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
		
		
		// measure total and roi signal across all slices
		for (slice=1; slice<(number_slices+1); slice++) {
			
			// measure roi signal
			roiManager("Select", 0);
			setSlice(slice);
			run("Measure");
			
			// store roi values in arrays
			roiMean_results[slice - 1] = getResult("Mean", 0);
			roiArea_results[slice - 1] = getResult("Area", 0);
			
			close("Results");
			
			// measure total signal
			roiManager("Deselect");
			selectImage("in_use");
			run("Measure");
			
			// store total values in arrays
			name_results[slice-1] = strip_name;
			slice_results[slice-1] = slice;
			time_results[slice-1] = (slice-1)*frame_interval;
			totalMean_results[slice-1] = getResult("Mean", 0);
			totalArea_results[slice-1] = getResult("Area", 0);
			close("Results");
		}
		
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
		
		// close everything
		close("saveTable");
		close("in_use");
		close("ROI Manager");
	}
}

// end
beep();