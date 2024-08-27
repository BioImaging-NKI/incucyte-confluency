// Macro to create a multiplot montage from analyzed Incucyte confluency data

//Layout of the plate
columns = 12;
rows = 6;

frameInterval = 4;
axisFontSize = 32;

table = "Averaged Confluency";		//The title of the input table
skipColumns = 2;					//skip 'Date Time' and 'Elapsed' columns in the input table

if(!isOpen("Averaged Confluency")) table = table+".tsv";
selectWindow(table);
timeArray = Array.getSequence(Table.size);
timeArray = multiplyArraywithScalar(timeArray, frameInterval);
//Plot.setXYLabels("Distance (microns)", "Gray Value");

setBatchMode(true);

plotString = "";
headings = Table.headings;
headers = split(headings, "\t");
for(i = skipColumns; i < headers.length; i++) {
	Plot.create(headers[i], "time (s)", "confluency (%)");
	Plot.setFrameSize(400, 400);
	Plot.setAxisLabelSize(axisFontSize);
	Plot.setFontSize(axisFontSize, "bold");
	Plot.setLineWidth(8);
	Plot.setLimits(0, 100, 0, 100);
	confluency = Table.getColumn(headers[i], table);
	Plot.add("connected", timeArray, confluency);
	Plot.setStyle(0, "blue,#7777ff,5,Connected Circles");
	Plot.setFontSize(axisFontSize*2, "bold");
	Plot.setFormatFlags("11000000000000");
	Plot.setXYLabels("", "");
	Plot.setColor("black");
	Plot.addText(headers[i], 0.05, 0.15);
	plotString += " image"+i+1-skipColumns+"=["+headers[i]+"]";
	Plot.show();
}
run("Concatenate...", "  title=All_plots open "+plotString);

//Crop the plots
makeRectangle(129, 29, 403, 403);		//crop all borders
//makeRectangle(120, 20, 422, 422);		//leave some border
run("Crop");

run("Make Montage...", "columns="+columns+" rows="+rows+" scale=0.5");
close("All_plots");
setBatchMode("show");


//Multiplies all elements of an array with a scalar and returns the new array
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}
