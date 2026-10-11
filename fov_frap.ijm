////// for measuring total fov signal across frap timepoints
//// open .tf2 file for each frap experiment
//// auto-adjust contrast on first frame
//// set lookup table to appropriate color
//// save as .tif

//// run macro
// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// set parameters
run("Set Measurements...", "area mean min redirect=None decimal=4");

// run macro in batch mode
setBatchMode(true);
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	
	// open whole fov image
	if (matches(path, ".*_\\d+\\.tif$")) {

		// open and auto-ajdust stack
		strip_name = replace(name, ".tif", "");
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
		fovMean_results = newArray(number_slices);
		fovArea_results = newArray(number_slices);
		
		// measure fov signal across all slices
		for (slice=1; slice<(number_slices+1); slice++) {
			setSlice(slice);
			run("Measure");
			
			// store total values in arrays
			name_results[slice-1] = strip_name;
			slice_results[slice-1] = slice;
			time_results[slice-1] = (slice-1)*frame_interval;
			fovMean_results[slice-1] = getResult("Mean", slice-1);
			fovArea_results[slice-1] = getResult("Area", slice-1);
		}
		close("Results");
		
		// store arrays in table and save
		Table.create("saveTable");
		Table.setColumn("file", name_results);
		Table.setColumn("slice", slice_results);
		Table.setColumn("time", time_results);
		Table.setColumn("fovMean", fovMean_results);
		Table.setColumn("fovArea", fovArea_results);
		Table.save(dir + "/" + strip_name + "_fov.csv");
		
		// close everything
		close("saveTable");
		close("in_use");
	}
}

// end
beep();