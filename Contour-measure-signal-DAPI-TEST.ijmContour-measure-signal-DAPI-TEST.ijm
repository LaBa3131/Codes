// TEST VERSION - Single file with step-by-step validation
// 1- Scale the file : scaling size of pixels in real life 
// 2- Substract the background 
// 3- Thresholding your image: tell the software wich pixels to quanitfy and which to ignore automatic , test Li, otsu, triangle black ,Hunag 
// 4- Define your measurments (set measurments)
// 5- To actually emasure the intensity 

// INSTRUCTIONS FOR TEST RUN:
// - Select one image from each of the three folders when prompted
// - The macro will pause at each step so you can inspect the results
// - Press OK in the dialog to continue to the next step
// - Close images manually if you want, or they'll be cleared at the end

// tidy up screen 
Tstart=getTime();
if (nImages>0) run("Close All"); // if there are 1 or more images - close all

if (isOpen("Results")) { // close results window, if open
	selectWindow("Results");
	close("Results");
}
print("\\Clear"); // empty log window
roiManager("reset"); // empty ROI manager
run("Options...", "iterations=1 count=1 black edm=32-bit"); // set Binary Options

//File opening 
setBatchMode(0); //0 makes images visible ; 1 hides them
openfolder1= getDirectory("Choose your channel of interest file");
openfolder2= getDirectory("Choose where DAPI is");
openfolder3= getDirectory("Choose where your outline is");
	ListArray1 = getFileList(openfolder1);
	ListArray2 = getFileList(openfolder2);
	ListArray3 = getFileList(openfolder3);
	print("Number of processed files",ListArray1.length,openfolder1);
	print("Number of processed files",ListArray2.length,openfolder2);
	print("Number of processed files",ListArray3.length,openfolder3);
	
	// Validate that all folders have the same number of files
	if (ListArray1.length != ListArray2.length || ListArray1.length != ListArray3.length) {
		exit("Error: Folders have different numbers of files!");
	}

// TEST MODE: Process only the first file
for (i = 0; i < 1; i++) {  // Changed from ListArray1.length to 1 for single file test
		print("\\n=== PROCESSING FILE " + (i+1) + " ===");
		
		// STEP 1: Open and process outline
		print("STEP 1: Opening outline image...");
		filename3=ListArray3[i];
		open(openfolder3+filename3);
		waitForUser("STEP 1: Outline image opened. Check the image and press OK to continue.");
		
		print("STEP 2: Subtracting background from outline...");
		run("Subtract Background...", "rolling=500");
		waitForUser("STEP 2: Background subtracted. Press OK to continue.");
		
		print("STEP 3: Thresholding outline...");
		setAutoThreshold("Triangle dark");
		waitForUser("STEP 3: Threshold applied. Press OK to continue.");
		
		print("STEP 4: Creating selection and adding to ROI Manager...");
		run("Create Selection");
		run("ROI Manager...");
		roiManager("Add");
		waitForUser("STEP 4: Outline ROI added. Press OK to continue.");
		
		// STEP 5: Open and process channel of interest
		print("STEP 5: Opening channel of interest image...");
		filename1=ListArray1[i];
		open(openfolder1+filename1);
		waitForUser("STEP 5: Channel image opened. Press OK to continue.");
		
		print("STEP 6: Converting to 16-bit and subtracting background...");
		run("16-bit");
		run("Subtract Background...", "rolling=1000");
		setOption("ScaleConversions", true);
		waitForUser("STEP 6: 16-bit conversion and background subtraction done. Press OK to continue.");
		
		print("STEP 7: Setting scale...");
		run("Set Scale...", "distance=1 known=0.325 unit=um");
		waitForUser("STEP 7: Scale set to 0.325 um/pixel. Press OK to continue.");
		
		print("STEP 8: Thresholding channel image...");
		setAutoThreshold("Triangle dark");
		waitForUser("STEP 8: Threshold applied. Press OK to continue.");
		
		print("STEP 9: Selecting outline ROI and measuring...");
		roiManager("Select", 0);
		run("Set Measurements...", "area mean integrated area_fraction limit display redirect=None decimal=2");
		run("Measure");
		waitForUser("STEP 9: Measurement complete. Check Results window. Press OK to continue.");
		
		// STEP 10: Open and process DAPI
		print("STEP 10: Opening DAPI image...");
		filename2=ListArray2[i];
		open(openfolder2+filename2);
		waitForUser("STEP 10: DAPI image opened. Press OK to continue.");
		
		print("STEP 11: Background subtraction and scaling...");
		run("Subtract Background...", "rolling=300");
		setOption("ScaleConversions", true);
		run("Set Scale...", "distance=1 known=0.325 unit=um");
		waitForUser("STEP 11: Background subtracted and scale set. Press OK to continue.");
		
		print("STEP 12: Thresholding DAPI...");
		setAutoThreshold("Triangle dark");
		setOption("BlackBackground", true);
		waitForUser("STEP 12: Threshold applied. Press OK to continue.");
		
		print("STEP 13: Converting to mask and watershed...");
		run("Convert to Mask");
		run("Watershed");
		waitForUser("STEP 13: Mask created and watershed applied. Press OK to continue.");
		
		print("STEP 14: Analyzing particles (nuclei detection)...");
		run("Analyze Particles...", "size=50-5000 circularity=0.10-1.00 summarize overlay add");
		waitForUser("STEP 14: Particle analysis complete. Check detected nuclei. Press OK to continue.");
		
		print("STEP 15: Recording cell number...");
		Nucnumber=roiManager("count") - 1; // Subtract the outline ROI
		setResult("cell number", nResults-1, Nucnumber);
		updateResults();
		print("Cell number recorded: " + Nucnumber);
		waitForUser("STEP 15: Cell number = " + Nucnumber + " Press OK to finish test.");
		
		// Clean up for next iteration (or end)
		close("*");
		roiManager("reset");
}

print("\\n=== TEST RUN COMPLETE ===");
print("All steps executed successfully!");
