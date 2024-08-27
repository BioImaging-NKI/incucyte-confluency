/* Macro to quantify the cell confluency in multiwell plates, imaged in the Incucyte system.
 * Input: two folders containing (paired) fluorescence and brightfield images, with exactly the same names.
 * Image names should ideally have the structure '...[wellName]_[position]_[yyyy]y[mm]m[dd]d_[hh]h[mm]m'.
 * 
 * Brief workflow:
 * - Fluorescence and brightfield images are loaded per well as image stacks.
 * - Background is subtracted from the fluorescence images (if applicable). The mean value of the background image is added again to prevent negative gray values. 
 *   This procedure proved more robust than a 'real' flatfield correction, for our images.
 * - Fluorescence image stacks are smoothed by a Gaussian blur
 * - Automatic threshold (Li) is calculated on the second frame (brightness of the 1st frame is slightly different),   
 *   using a histogram that is clipped at 'upperDisplaySetting' (the no-reset option), for more robustness. This threshold is then applied to the stack.
 *   A bias from this automatic threshold can be set in the script parameter dialog.
 * - A user-set number of Opening and Dilation steps are applied to the mask.
 * - Area percentages are measured and saved into a table (Confluency all positions).
 * - Results from different positions in each well are averaged and saved in a new table (Averaged confluency). 
 * - Composite image stacks from each position can be saved (binned). Channels: fluorescence-Mask-brightfield
 *  
 * Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 */

#@ File (label = "Input folder fluorescence images", style = "directory", description="Files in both folders should be named the same.") inputFluoFolder
#@ File (label = "Input folder brightfield images", style = "directory", description="Files in both folders should be named the same.") inputBrightfieldFolder
#@ File (label = "Output folder", style = "directory") outputFolder
#@ Boolean (label = "Correct background?", value = true) subtractBackground
#@ File (label = "...with background image", style = "file", description="If present, fluorescence images will be background-corrected", value="N/A") backgroundImagePath
#@ Integer (label = "Upper display setting for thresholding (~intensity of bright cells)", min=1, value=20) upperDisplaySetting
#@ Double (label = "Time interval between acquitisions (h)", min=0, value=4) frameInterval
#@ Double (label = "Sigma of Gaussian blur before thresholding", min=0, value=3) GaussianBlurSigma
#@ Double (label = "Threshold bias (shift from Auto threshold)", style="format:0.0", value=0.5) thresholdBias
#@ Integer (label = "nr. of opening iterations", min=0, value=2) openingIterations
#@ Integer (label = "nr. of dilation iterations", min=0, value=2) dilationIterations
#@ String (label = "String ending at the well name (e.g. \"well_\")", value = "well_") stringBeforeWellName
#@ Boolean (label = "Save composite images (fluorescence, mask, brightfield)", value=true) saveImages
#@ Integer (label = "To save disk space, xy-bin images with factor", min=1, value=2) binFactor

version = 1.0;

useDateAndTime = true;	//Image names should have the structure '..._2022y05m05d_16h46m'. Set to false if the macro crashes because of this.
multiplication = 10;	//Scale up factor for the fluorescence channel in the saved composite images (saved as 8-bit, but need some multiplication to preserve dynamic range).
maxDisplaySetting = 120;//upper display setting in the saved images.


setOption("ExpandableArrays", true);
run("Conversions...", " ");
run("Set Measurements...", "area area_fraction limit redirect=None decimal=3");

//Create parameter list
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
List.set("00. ==Confluency Analyzer settings=", "");
List.set("01. Date ", " " + DayNames[dayOfWeek] + " " + dayOfMonth + " " + MonthNames[month] + " " + year);
List.set("02. Time ", " " + hour +"h"+IJ.pad(minute,2));
List.set("03. Version ", version);
List.set("04. Input fluorescence folder ", inputFluoFolder);
List.set("05. Input brightfield folder ", inputBrightfieldFolder);
List.set("06. Output folder ", inputFluoFolder);
List.set("07. Correct background? ", subtractBackground);
List.set("08. With background image ", backgroundImagePath);
List.set("09. Upper display setting for thresholding (intensity of bright cells) ", upperDisplaySetting);
List.set("10. Time interval between acquitisions (h) ", frameInterval);
List.set("11. Sigma of Gaussian blur before thresholding ", GaussianBlurSigma);
List.set("12. Threshold bias (shift from Auto threshold) ", thresholdBias);
List.set("13. nr. of opening iterations ", openingIterations);
List.set("14. nr. of dilation iterations ", dilationIterations);
List.set("15. String ending at the well name (e.g. \"well_\") ", stringBeforeWellName);
List.set("16. Save composite images (fluorescence, mask, brightfield) ", dilationIterations);
List.set("17. To save disk space, xy-bin images with factor ", stringBeforeWellName);
list = List.getList();
logFile = File.open(outputFolder + File.separator + "Analysis settings -"+List.getValue("01. Date ") + List.getValue("02. Time ")+".txt");
print(logFile, list);

if(!File.exists(outputFolder)) File.makeDirectory(outputFolder);
fileList = getFileList(inputFluoFolder);
indexOfWell = lastIndexOf(fileList[0], stringBeforeWellName)+stringBeforeWellName.length;

print("\\Clear");

//Get existing well names and positions
showStatus("Finding all well names...");
wells = newArray();
for(i=0; i<fileList.length; i++) {
	showProgress(i, fileList.length);
	wells[i] = substring(fileList[i], indexOfWell, indexOfWell+1) + IJ.pad(substring(fileList[i], indexOfWell+1, indexOf(fileList[i], "_", indexOfWell)),2);
}
uniqueWellNames = uniqueValuesInArray(wells);

//Sort names, then de-pad.
uniqueWellNames = Array.sort(uniqueWellNames);
for(i=0; i<uniqueWellNames.length; i++) {
	if(substring(uniqueWellNames[i], 1, 2) == "0") uniqueWellNames[i] = substring(uniqueWellNames[i], 0, 1) + substring(uniqueWellNames[i], 2, uniqueWellNames[i].length);
}

positions = newArray();
for(i=0; i<fileList.length; i++) {
	positions[i] = substring(fileList[i], indexOf(fileList[i], "_", indexOfWell)+1, indexOf(fileList[i], "_", indexOf(fileList[i], "_", indexOfWell)+1));
}
uniquePositions = uniqueValuesInArray(positions);

nrWells = uniqueWellNames.length;
nrPos = uniquePositions.length;

print("Found "+nrWells+" wells with "+nrPos+" positions.");

setBatchMode(true);
Table.create("Thresholds");
ct = "Confluency all positions";
Table.create(ct);
firstRun = true;
for(w=0; w<nrWells; w++) {
	showProgress(w, nrWells);
	for(p=0; p<nrPos; p++) {
		run("Close All");
		wellname = uniqueWellNames[w]+"_"+uniquePositions[p];
		print("Processing "+wellname+"...");
		showStatus("Processing "+wellname+"...");
		run("Clear Results");
		File.openSequence(inputFluoFolder, " filter=_"+wellname);
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		image = getTitle();
		run("32-bit");
		run("Grays");
		getDimensions(width, height, channels, slices, frames);

		if(firstRun == true && useDateAndTime == true) dateAndTimeArray = getTimeStamp(image);
		
		//Subtract background (background image has to be previously constructed) - not a real flatfield correction, but a subtraction!
		if(subtractBackground == true) {
			if(File.exists(backgroundImagePath)) open(backgroundImagePath);
			else exit("Please provide a valid path to a background image.");
			background = getTitle();
			meanBackground = getValue("Mean");
			imageCalculator("subtract 32-bit stack", image, background);
			run("Add...", "value="+meanBackground+" stack");	//Add the average of what was subtracted, to prevent values below zero
			selectWindow(image);
		}
		run("Enhance Contrast", "saturated=0.35");
		run("Duplicate...", "title=Mask duplicate");

		run("Gaussian Blur...", "sigma="+GaussianBlurSigma+" stack");
		if(frames>1) Stack.setFrame(2);	//frame 1 is somehow more bright. Unclear whether this is always the case, but better safe than sorry.

		setMinAndMax(0, upperDisplaySetting);
		setAutoThreshold("Li dark no-reset");
		getThreshold(lower, upper);
		setThreshold(lower + thresholdBias, upper);
		index = w*nrPos + p;
		Table.set("well_pos", index, wellname, "Thresholds");
		Table.set("threshold", index, lower + thresholdBias, "Thresholds");
		Table.update;
		run("Convert to Mask", "method=Li background=Dark black");

		run("Options...", "iterations="+openingIterations+" count=1 black do=Nothing");
		run("Open", "stack");
		run("Options...", "iterations="+dilationIterations+" count=1 black do=Nothing");
		run("Dilate", "stack");

		setThreshold(1, 255);
		run("Measure Stack...");
		areaFraction = Table.getColumn("%Area", "Results");
		Table.setColumn(wellname, areaFraction, ct);

		//Multiply fluorescence intensity image with [multiplication] and scale to 8-bit - Modify the numbers at the top of the script if your images have a very different range 
		selectWindow(image);
		run("Multiply...", "value="+multiplication);
		setMinAndMax(0, 255);
		run("8-bit");

		//Open brightfield images, merge with fluorescence images and the calculated mask
		if(saveImages == true) {
			File.openSequence(inputBrightfieldFolder, " filter=_"+wellname);
			brightfield = getTitle;
			if(bitDepth() != 8) run("8-bit");
			run("Merge Channels...", "c1="+image+" c2=Mask c3="+brightfield+" create");
			Stack.setChannel(3);
			run("Grays");
			run("Enhance Contrast", "saturated=0.35");
			Stack.setChannel(2);
			run("Blue");
			Stack.setChannel(1);
			run("Red");
			Stack.setFrame(1);
			median = getValue("Median");
			setMinAndMax(median, maxDisplaySetting);
			run("Bin...", "x="+binFactor+" y="+binFactor+" z=1 bin=Average");
			setMetadata("info", list);
			saveAs("Tiff", outputFolder + File.separator + wellname);
		}
		close();
	}
}

selectWindow("Thresholds");
Table.save(outputFolder + File.separator + "Thresholds -"+List.getValue("01. Date ") + List.getValue("02. Time ")+".tsv");

selectWindow(ct);
Table.save(outputFolder + File.separator + ct + " - " + List.getValue("01. Date ") + List.getValue("02. Time ")+".tsv");

//Cumbersome and slow method of averaging Table columns, but it works.
//A much fast alternative would be to create an image, bin in x-direction, and get the data.
avgct = "Averaged Confluency";
Table.create(avgct);
elapsed = Array.getSequence(frames);
elapsed = multiplyArraywithScalar(elapsed, frameInterval);
if(useDateAndTime) Table.setColumn("Date Time", dateAndTimeArray, avgct);
Table.setColumn("Elapsed", elapsed, avgct);

headings = Table.headings(ct);
headers = split(headings, "\t");
for(i = 0; i < headers.length; i=i+nrPos) {
	names = Array.slice(headers, i, i+nrPos);
	avgConfluency = averageTableColumns(names, ct);
	Table.setColumn(substring(headers[i], 0, indexOf(headers[i], "_")), avgConfluency, avgct);
}
Table.update;
Table.save(outputFolder + File.separator + avgct + " - "+List.getValue("01. Date ") + List.getValue("02. Time ")+".tsv");


//Retrieves the date and time from the frame labels as an array of strings
function getTimeStamp(image) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	dateAndTimeArray = newArray(frames);
	
	end = newArray(frames);
	min = newArray(frames);
	hour = newArray(frames);
	day = newArray(frames);
	month = newArray(frames);
	year = newArray(frames);

	for(f=1; f<=frames; f++) {
		Stack.setFrame(f);
		string = getMetadata("label");
		end = lastIndexOf(string, ".tif");
		minute = substring(string, end-3, end-1);
		hour = substring(string, end-6, end-4);
		day = substring(string, end-10, end-8);
		month = substring(string, end-13, end-11);
		year = substring(string, end-18, end-14);
		
		dateAndTimeArray[f-1] = day + "-" + month + "-" + year + " " + hour + ":" + minute;
	}
	return dateAndTimeArray;	
}


//Multiplies all elements of an array with a scalar and returns the new array
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}


//Returns true if the value occurs within the array
function occursInArray(array, value) {
	for(i=0; i<array.length; i++) {
		if(array[i] == value) return true;
	}
	return false;
}


//Returns an array with unique values in the array
function uniqueValuesInArray(array) {
	count=0;
	uniqueArray = newArray(array.length);
	for (a=0; a<lengthOf(array); a++) {
		if(!occursInArray(Array.trim(uniqueArray, count), array[a]) && !matches(array[a],"None")) {
			uniqueArray[count]=array[a];
			count++;
		}
	}
	return Array.trim(uniqueArray, count);
}


//Average multiple columns of a table element-wise. 'array' is a string array containing the array variable names
function averageTableColumns(nameArray, table) {
	selectWindow(table);
	summed_Array = newArray(Table.size);
	for(i=0; i<Table.size; i++) {
		sum = 0;
		for (k = 0; k<nameArray.length; k++) {
			sum += Table.get(nameArray[k], i);
		}
		summed_Array[i] = sum / nameArray.length;
	}
	return summed_Array;
}
