//Macro to create a background file from many individual images, by making a median projection.
//The idea behind it is that most of the time a pixel is npt occupied by a cell. Works best for images/wells that are not fully overgrown.
 
#@ File (label = "Input folder containing single plane images", style = "directory") inputFolder
#@ Integer (label = "Number of images to load") nrImages

File.openSequence(inputFolder, "count="+nrImages);
run("Z Project...", "projection=Median");
dir = File.getParent(inputFolder);
saveAs("Tiff", dir + File.separator + "background");
