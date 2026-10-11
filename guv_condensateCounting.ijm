//// don't denoise for frap analysis but consider it for representative images
//// open .tf2 file for each frap experiment
//// auto-adjust contrast on first frame
//// set lookup table to appropriate color
//// save as .tif
//// duplicate and crop guv of interest (rectangle that's twice as large as guv)
//// auto-adjust contrast
//// save as "_frap#.tif"
/// draw segmented line roi (30 pixels thick) on guv membrane
/// save roi as "_frap#_circle.roi"

//// run macro
// select directory
dir = getDirectory("Choose Source Directory ");
list = getFileList(dir);

// set measurements
run("Set Measurements...", "area mean min redirect=None decimal=4");

// run macro in batch mode
setBatchMode(true);
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = dir + list[i];
	name = File.getName(path);
	
	// process circle roi
	if (endsWith(path, "_circle.roi")) {

		// open and auto-ajdust stack
		strip_name = replace(name, "_circle.roi", "");
		open(strip_name + ".tif");
		rename("in_use");
		setSlice(1);
		resetMinAndMax();
		run("Enhance Contrast", "saturated=0.35");
		
		// open and straighten roi
		open(strip_name + "_circle.roi");
		run("Straighten...", "title=straight line=30");
		selectImage("straight");
		getPixelSize(unit, pixelWidth, pixelHeight);
		
		// caculate box dimensions
		numCols = getWidth();
		numRows = 1;
		boxWidth = 1;
		boxHeight = getHeight();
		
		// make arrays
		name_results = newArray(numCols);
		boxNumber_results = newArray(numCols);
		boxRightPixel_results = newArray(numCols);
		boxRightMicron_results = newArray(numCols);
		boxHeightPixel_results = newArray(numCols);
		boxHeightMicron_results = newArray(numCols);
		boxAreaPixel_results = newArray(numCols);
		boxAreaMicron_results = newArray(numCols);
		boxMean_results = newArray(numCols);
		
		// cycle through guv membrane boxes to quantify intensities
		for (row = 0; row < numRows; row++) {
  			for (col = 0; col < numCols; col++) {
    			
 		     	// calculate box coordinates
 		     	boxX = col;
  		    	boxY = row;

   		   		// define ROI for current box
   		   		makeRectangle(boxX, boxY, boxWidth, boxHeight);

    		  	// calculate box intensity
    		  	run("Measure");
			}			
		}
		
		// store values in arrays
		for(col = 0; col < numCols; col++) {
			name_results[col] = strip_name;
			boxNumber_results[col] = col;
			boxRightPixel_results[col] = col + 1;
			boxRightMicron_results[col] = (col + 1)*pixelWidth;
			boxHeightPixel_results[col] = boxHeight;
			boxHeightMicron_results[col] = boxHeight*pixelHeight;
			boxAreaPixel_results[col] = boxWidth*boxHeight;
			boxAreaMicron_results[col] = boxWidth*pixelWidth*boxHeight*pixelHeight;
			boxMean_results[col] = getResult("Mean", col);
		}
		close("Results");
		
		// store arrays in table and save
		Table.create("saveTable");
		Table.setColumn("file", name_results);
		Table.setColumn("boxNumber", boxNumber_results);
		Table.setColumn("boxRightPixel", boxRightPixel_results);
		Table.setColumn("boxRightMicron", boxRightMicron_results);
		Table.setColumn("boxHeightPixel", boxHeightPixel_results);
		Table.setColumn("boxHeightMicron", boxHeightMicron_results);
		Table.setColumn("boxAreaPixel", boxAreaPixel_results);
		Table.setColumn("boxAreaMicron", boxAreaMicron_results);
		Table.setColumn("boxMean", boxMean_results);
		Table.save(dir + "/" + strip_name + "_circle.csv");
		
		// close everything
		close("saveTable");
		close("in_use");
		close("straight");	
	}
}

// end
beep();