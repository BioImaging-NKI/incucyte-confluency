# Incucyte multiwell plate confluency analysis
ImageJ macros for the analysis of fluorescence timelapse multiwell plate images from the Incucyte system, developed for the paper "...", submitted to "...".

In order to use this script, fluorescence and brightfield timelapses from the Incucyte system should be exported as single-plane `tif` images with wellname, position and timestamp in the file name. (an example file name is `well_B1_2_2022y05m08d_08h46m.tif`). Fluorescence and brightfield images should have exactly the same name, in separate folders.  The ImageJ1 macro `Incucyte confluency analysis.ijm` combines images for each position, preprocesses them and analyses the confluency over time. Measurements are saved in a table `Confluency all positions`, as well as in a table `Averaged Confluency`, where all positions in a well have been averaged.

## Installation
Download the three `.ijm` files. Drag&drop a file onto the Fiji window opens the script editor. Click `Run` to run the macro. :-)

## Brief workflow
- Fluorescence and brightfield images are loaded per well as image stacks. The brihgtfield images are not used in the quantifications, 
- Background can be subtracted from the fluorescence images. If you don't have an background image, or 'empty' image, the macro `Incucyte confluency analysis - create background.ijm` can be used to generate one. The mean value of the background image is added again to prevent negative gray values. This procedure proved more robust than a 'real' flatfield correction, for our low-intensity images.
- Fluorescence image stacks are smoothed by a Gaussian blur with user-defined sigma.
- Automatic threshold (Li) is calculated on the second frame (the brightness of the 1st frame was found to be slightly different), using a histogram that is clipped at 'upperDisplaySetting' (the no-reset option), for more robustness. This threshold is then applied to the stack.
- A (absolute) bias from this automatic threshold can be set in the script parameter dialog to tweak the automatic threshold
- A user-set number of Opening and Dilation steps are applied to the masks, to remove single pixels and expand the masks.
- Area percentages of the masks over time are measured and saved into a table (Confluency all positions).
- Results from different positions in each well are averaged and saved in a new table (Averaged confluency).
- Composite image stacks from each position can be saved for visual inspection (with user-defined binning). Channels: ch1=fluorescence (red), ch2=Confluency mask (blue), ch3=brightfield (grays).

## Output examples

A background image can be generated by taking for each pixel position the median for all images in a folder. This works well with many images with relatively low confluency. Subtracting the background will help avoiding artefacts in the thresholding and ultimately in the measured confluency, especially when the fluorescence signal is weak.
![image](https://github.com/user-attachments/assets/a6324dc3-a87a-47fd-b387-3ba5679f64c9)

In the paper cells were co-cultured with (unstained) cytotoxic T cells. Below the resulting images of two different conditions. In the right image cells are killed by the T cells before they can grow out.
![confluency_2wells](https://github.com/user-attachments/assets/4c2d107f-82ab-454b-8a64-f5e2aa4761d6)

After running an analysis, the resulting `Average Confluency` table should look like this:
![image](https://github.com/user-attachments/assets/841debb9-f6c4-41cf-94c7-6ef4b077fa04)

Finally, the macro `Incucyte confluency analysis - create multiplot from table` can be used to generate an overview of confluency over time for the whole well plate. It reads the Average Confluency table and produces a montage image:
![image](https://github.com/user-attachments/assets/16478fef-e6f4-43e9-8441-33ae8971aa01)

