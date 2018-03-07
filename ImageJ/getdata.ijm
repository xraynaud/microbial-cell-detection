setBatchMode(true)

//basedir = "/home/xraynaud/Mycore/Projets/Rhizoplane_Dist/Data_for_paper/Haze experiment/";


basedir =  getDirectory("Choose root directory (Should contain ./Images and ./Data)"); //"/home/xraynaud/Mycore/Projets/Rhizoplane_Dist/Data/Gnoto_Alex/new/";

list = getFileList(basedir+"Images/");

for (i = 0; i <list.length; i++) {
	if (endsWith(list[i],".tif")) {
		print(list[i]);
		open(basedir+"Images/"+list[i]);
		analyse_image(true,0,0,1,basedir+"Data/",false);
		close(list[i]);
	}
}

function analyse_image(props,FOV,walls,cells,dir, verbose) {

// This function takes a multichannel composite stack and makes the data extraction as described in
// Schmidt et al. 
// This function should be copied in file Xxx to be available to other ImageJ/Fiji macro.
// Function parameters are:
// 	- props (Boolean): create a R source file (filename_props.R) containing various properties of the image
//		(size, resolution...)
// 	- FOV (numeric): image channel to calculate a ROI in which low signal parts of the image has been 
//		removed.Produce a text file (filename_fov.txt) that can be imported into R as a spatstat owin() 
//		object. 0 will use the full image size as ROI.
//	- walls (numeric): do the cell walls extraction on given channel. Requires the Ridge Detection plugin 
//		(https://github.com/thorstenwagner/ij-ridgedetection) to work.  Produces a TIFF file 
//		(filename_segments.tiff) that can be converted to SVG using autotrace 
//		(http://autotrace.sourceforge.net/) or any other vectorizing software that can do centerline tracing.
//		0 does not do wall extraction.
// 	- cells (numeric): do the cell coordinate extraction on given channel. Requires the Ridge Detection 
//		plugin (https://github.com/thorstenwagner/ij-ridgedetection) to work. Produces a _raw_ coordinate 
//		set that must be cleaned using R function XXX. 0 does not do cell extraction.
// 	- dir (text): full path of directory where to save files. Error if path does not exist
// 	- verbose (Boolean): should intermediate images as well as text output be shown.

//make macro operate in background if verbose false.
	//setBatchMode(!verbose);

// get some info on the image
	img=getTitle; 
	name = replace(img,".tif","");
	getDimensions(w, h, nchannels, slices, frames);	     
 	getVoxelSize(size_x, size_y, size_z, unit);
	if (verbose) {
		print(name + " : " +slices+" slices in "+ nchannels +" channels");
	}
	path_fov = dir+"/"+name+"_fov.txt";
// reset ROI manager and overlay.
	run("Remove Overlay");
 	roiManager("reset");
 	run("Clear Results");

// Separate image channels, make a maximum intensity projection for each channel. New images are named ch1, ch2...
	for (c=0; c<nchannels;c++) {
		selectWindow(img);
		run("Duplicate...", "title=ch"+(c+1)+" duplicate channels="+(c+1));
  		selectWindow("ch"+(c+1));
		run("Z Project...", "projection=[Max Intensity]");
	}

//Image Properties
	if (props) {
		path_props = dir+"/"+name+"_props";
		
		props_file = File.open(path_props+".txt");
		print(props_file,"array = c("+w+","+h+","+slices+")");
		print(props_file,"xrange = c(0,"+w*size_x+")");
		print(props_file,"yrange = c(0,"+h*size_y+")");
		print(props_file,"zrange = c(0,"+slices*size_z+")");
		File.close(props_file);
		if (File.rename(path_props+".txt",path_props+".R")) {
			if (verbose){
				print("Image properties saved to file "+path_props+".R");
			}
		} else {
			print("Cannot save file propoerties to "+path_props+".R");
		}
		selectWindow(img);
		run("Z Project...", "projection=[Max Intensity]");
		run("RGB Color");
		path_zproj = name+"_maxint.tiff";
		saveAs("TIF", dir+"/"+path_zproj);
		close("MAX_"+img);
		close(path_zproj);
		if (verbose) {
			print("Saved a RGB z-projection image to "+dir+"/"+path_zproj);
		}
	}

//FOV
	fov = File.open(path_fov);
	print(fov,"id\tx\ty");
	if (FOV>0) {
		selectWindow("MAX_ch"+FOV);
		if (verbose) {
			print("Analysing Field of view from channel "+FOV);
		}
		run("Duplicate...", "title=FOV");
		run("Subtract Background...", "rolling=10 create disable");
		selectWindow("FOV");
		run("Median...", "radius=5");
		//run("Gaussian Blur...", "sigma=5");
		getHistogram(0,counts,256);
		minbin =0;
		for (i=0;i<256;i++) {
			if (minbin == 0 && counts[i] > 0) {
				minbin = i;
			}
		}
		setThreshold(minbin+5, 255);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Fill Holes");
		run("Analyze Particles...", "size=50-Infinity show=Outlines add in_situ");
		lengths = newArray(roiManager("count"));
		regions = roiManager("count");
		
		for (i=0; i<regions; i++) {
			roiManager("select", i);
			getSelectionCoordinates(x_, y_); 
			x_ = Array.concat(x_,x_[0]);
			y_ = Array.concat(y_,y_[0]);
			toScaled(x_,y_);
			for (p=0; p < x_.length;p++) {
				print(fov,(i+1)+"\t"+x_[p]+"\t"+y_[p]);
			}
		}
		close("FOV");
		if (!verbose) {
			roiManager("reset");
		}
	} else {
		if (verbose) {
			print("Field of view is image size");
		}
		print(fov,"1\t0\t0");
		print(fov,"2\t0\t"+h*size_y);
		print(fov,"3\t"+w*size_x+"\t"+h*size_y);
		print(fov,"4\t"+w*size_x+"\t0");
	}
	File.close(fov);
	if (verbose) {
		print("ROI of field of view saved to "+ path_fov);
	}

//Walls
	if (walls >=0) {
		if (walls == 0) {
			chmax = nchannels ;
			chmin = 1; 
			if (verbose) {
				print("Extracting root cell walls from all channels");
			} 
		} else {
			chmax = walls;
			chmin = walls;
			if (verbose) {
				print("Extracting root cell walls from channel "+ chmax);
			}
		}
		
		
		for (c=chmin; c<=chmax;c++) {
			if (verbose) {
				print("analysing channel "+(c));
			}
  			selectWindow("ch"+(c));
  			run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
  			run("FeatureJ Laplacian", "compute smoothing=10");
  			run("Unsharp Mask...", "radius=3 mask=0.90 stack");
  			run("8-bit");
  			run("Ridge Detection", "line_width=15 high_contrast=200 low_contrast=0 darkline process_stack make_binary method_for_overlap_resolution=NONE sigma=4.83 lower_threshold=0 upper_threshold=0.68 minimum_line_length=15 maximum=0");
  			close("ch"+(c)+" Laplacian");
  			if (isOpen("ch"+(c)+" Laplacian Detected segments")) {
  				run("Z Project...", "projection=[Max Intensity]");
  			}
 		}
 		if (isOpen("MAX_ch1 Laplacian Detected segments") && isOpen("MAX_ch2 Laplacian Detected segments")) {
 			imageCalculator("Add create", "MAX_ch1 Laplacian Detected segments","MAX_ch2 Laplacian Detected segments");
 		}
 		run("Close-");
 		run("Skeletonize");
		path_segments = name+"_segments.tiff";
		saveAs("TIF", dir + "/" + path_segments);
		if (isOpen("Result of MAX_ch1 Laplacian Detected segments")) {
			close("Result of MAX_ch1 Laplacian Detected segments");
		}
		if (isOpen("MAX_ch1 Laplacian Detected segments")) {
 			close("MAX_ch1 Laplacian Detected segments");
 		}
 		if (isOpen("MAX_ch2 Laplacian Detected segments")) {
 			close("MAX_ch2 Laplacian Detected segments");
 		}
		//close("Noise");
		//close("Result of MAX_ch"+walls);
		//close("Mask of Noise");
		if (verbose) {
			print("Detected segments saved as "+ dir + "/" + path_segments);
		} else {
			close(path_segments);
		}
	}

// Cells
	if (cells) {
		if (verbose) {
			print("Analysing cell coordinates from channel " + cells);
		}
  		pts=0;
		for (s=1; s<=slices; s=s+1) { 
			selectWindow("ch"+cells);
			if (verbose) {
				print("Working on slice "+ s);
			}
			setSlice(s);
			run("Duplicate...","title=Slice use");
			selectWindow("Slice");
			getStatistics(area, mean, min, max, std, histogram);
			if (verbose) {
				print("Max signal in slice "+s+": "+max);
			}
		
			if (max > 0) { // Only work on slices that have signal
				if (walls>0) {
					selectWindow("ch"+walls);
					setSlice(s);
					run("Duplicate...","title=Other use");
				//	run("Median...", "radius=3");
					run("Subtract Background...", "rolling=10");
					imageCalculator("Subtract create", "Slice","Other");
					close("Other");
				} else {
					run("Duplicate...","title=[Result of Slice] use");
				}
				selectWindow("Result of Slice");
				//run("Subtract Background...", "rolling=10");
				run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");

				run("Duplicate...", "title=Slice_Threshold");
				selectWindow("Slice_Threshold");
				setMinAndMax(0, 80);
				run("Apply LUT");
	
				setAutoThreshold("Otsu dark");
				run("Convert to Mask");	
				//run("Median...", "radius=2");  // Cette ligne modifie subtilement les résultats. 

				if (walls >0) {
					if (verbose) {
						print("Looking for linear features");
					}
					selectWindow("Result of Slice");
					run("Ridge Detection", "line_width=6 high_contrast=40 low_contrast=0 make_binary method_for_overlap_resolution=NONE sigma=2.23 lower_threshold=0 upper_threshold=0.51 minimum_line_length=40 maximum=0");
					run("Dilate");
					run("Dilate");	
					imageCalculator("Subtract create", "Slice_Threshold","Result of Slice Detected segments");
					if (verbose) {
						print("Finished looking for linear features");
					}
				}
				run("Duplicate...","title=[Result of Slice_Threshold] use");
				selectWindow("Result of Slice_Threshold");
				//run("Median...", "radius=1");  // Cette ligne modifie subtilement les résultats. 

				run("Create Selection");
				selectType = selectionType();
				if (selectType>=0) {
					roiManager("Add");
					selectWindow("Slice");
					roiManager("Select", 0);	

					run("Find Maxima...", "noise=25 output=[Point Selection]");
		
					getSelectionCoordinates(x, y);
		//if (x.length>0) { 
					truex = Array.copy(x);
					truey = Array.copy(y);
					toScaled(truex,truey); 
	// truex, truey contient les coordonées des points contenus dans la couche actuell s

		   	 		for (i=0; i<x.length; i++) {
    		 			setResult("X",pts,truex[i]);
     					setResult("Y",pts,truey[i]);
     					setResult("Z",pts,(s-1)*size_z);
     					setResult("I", pts, getPixel(x[i],y[i]));
    		 			pts=pts+1;
   	 				}
    				if (verbose) {
    					print("Found "+x.length+" cells");
    				}
   		 			roiManager("Delete");
    				roiManager("Add");
    				if (selectionType() == 10) {
    					selectWindow(img);
						run("From ROI Manager");
						Overlay.setPosition(cells, s, 1)
    				}
    				roiManager("Delete");
				}
				saveAs("Results", dir+"/"+name+"_coords.csv");
				close("Results");
				close("Result of Slice");
				close("Result of Slice_Threshold");
				close("Result of Slice Detected segments");
				close("Slice_Threshold");
				close("Result of Slice_Threshold");
				close("Slice Detected segments");
				close("Slice");
				roiManager("reset");
			}
		}
	}
	for (c = 0; c < nchannels; c++) {
		close("ch"+(c+1));
		close("MAX_ch"+(c+1));
	}
	if (verbose) {
		print("End of analysis");
	}
// reset setbatchmode to normal.
	//setBatchMode(false);
}

	