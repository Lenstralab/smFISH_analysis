/*
Version:  See below
Description: This macro has the option to:
   * Open montage images created with Python's FISH pipeline
   * Select images and/or regions to delete in the cell and nucleus segmentation images in that directory

Contact:   Linda Joosen
Internet:  -
To do:
   * 
   * 
 */
	
//  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  
	title		=	"Montage_Delete";
	version		=	"v1";
	MACRO 		= 	title+version;
	date		=	"30 Jan 2020";
	Contact		= 	"linda_joosen@hotmail.com"
	Internet	= 	"https://spark.adobe.com/page/5jtPRBEYfeNjp/";
//  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  

macro MACRO {

	//	1  Start  ========================================================================================
	//	1.1  Start - Variables  ========================================================================================
	DeleteImages = "Yes";					// default for deleting images
	DeleteRegions = "Yes";					// default for deleting regions
	delimiter = ".tiff";					// to delete .tiff in filename, if present
	delimiter2 = "_montagePython_max.tif"; 	// to delete max montage name from in filename, for saving log file

	width = 2048;							// width of original images
	height = width;							// height of original images
	Scale = 4;								// scale of original- to montage image
	
	//	1.2  Start - Close (images, log, ROI)  ========================================================================================
	run("Close All");
	print("+~=-=~+");
	selectWindow("Log");
	print("\\Clear");
	roi = roiManager("count") - 1;
	if (roi >= 0) {
		roiManager("Deselect");
		roiManager("Delete");
	}

	//	1.3  Start - Dialog  ========================================================================================
	Dialog.create("User Dialog : Segmentation Image");    
		items = newArray("Outlined only", "Colored only", "B&W only", "Not", "Outlined & Colored", "Outlined & B&W");
	Dialog.addRadioButtonGroup("              Segmentation style viewed : ", items, 2, 3, "Outlined & Colored");
		items = newArray("Yes", "No");
	Dialog.addRadioButtonGroup("Show max projected channels (only when outlined only) : ", items, 1, 2, "Yes");
	Dialog.show();    
	LUT = Dialog.getRadioButton();
	MaxChan = Dialog.getRadioButton();

	StartTime = getTime();
	
	//	1.4  Start - Date & time  ========================================================================================
	print ("~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~"); 
	print("Macro :   ", title + ": " + version); 
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
		if (DayNames[dayOfWeek] == "Sun") { DAY = "Sunday"; }	if (DayNames[dayOfWeek] == "Mon") { DAY = "Monday"; }	if (DayNames[dayOfWeek] == "Tue") { DAY = "Tuesday"; }	if (DayNames[dayOfWeek] == "Wed") { DAY = "Wednesday"; }
		if (DayNames[dayOfWeek] == "Thu") { DAY = "Thursday"; }	if (DayNames[dayOfWeek] == "Fri") { DAY = "Friday"; }	if (DayNames[dayOfWeek] == "Sat") { DAY = "Saturday"; }
		TimeString ="Date: "+DayNames[dayOfWeek]+" ";
		if (dayOfMonth<10) {TimeString = TimeString+"0";}
		TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"     Time: ";
		if (hour<10) {TimeString = TimeString+"0";}
		TimeString = TimeString+hour+":";
		if (minute<10) {TimeString = TimeString+"0";}
		TimeString = TimeString+minute+":";
		if (second<10) {TimeString = TimeString+"0";}
		TimeString = TimeString+second;	
	print(TimeString);
		TimeStart = hour * 60;		TimeStart = TimeStart + minute * 60;	TimeStart = TimeStart + second;			
	// Date scientific
		DateString ="";
		if (month<10) {DateString = DateString+"0";}		Month = month + 1;						DateString = DateString+Month;
		if (dayOfMonth<10) {DateString = DateString+"0";}	DateString = DateString+dayOfMonth;
		Year = year - 2000;									YearString = "20" + Year;				DateString = YearString + DateString;
		//	print(DateString);	
	print ("Dimensions desktop screen : ", screenWidth, " * ", screenHeight);
	print ("Contact : ",Contact ,",   for more information  (check link in next line)");
	print(Internet);
	print ("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print ("* Files in directory : ");
	
	//	1.5  Start - Open directory  ========================================================================================
	filedir = getDirectory("Choose Source Directory ");

	//	1.6  Start - Define LUT positions & channels to add  ========================================================================================
	LUToutline = "No"; 	if (LUT == "Outlined only") { LUToutline = "Yes"; } 	else { if (LUT == "Outlined & Colored") { LUToutline = "Yes"; }  	else { if (LUT == "Outlined & B&W") { LUToutline = "Yes"; }}}
	LUTcolor = "No"; 	if (LUT == "Colored only") { LUTcolor = "Yes"; } 		else { if (LUT == "Outlined & Colored") { LUTcolor = "Yes"; }}
	LUTbw = "No"; 		if (LUT == "B&W only") { LUTbw = "Yes"; } 				else { if (LUT == "Outlined & B&W") { LUTbw = "Yes"; }}
	LUTbwo = "No"; 		if (LUTbw == "Yes") { 	if (LUToutline == "Yes") { LUTbwo = "Yes"; } }
	else { LUTbwo = "No"; 		if (LUTcolor == "Yes") { 	if (LUToutline == "Yes") { LUTbwo = "Yes"; } } }
	
	// Add channels
	ChannelAdd = 0;	// sets to 0
	if (LUTbwo == "Yes") {  // if both 2 segment type is chosen (e.g. color + contour  or  B&W + Contour)
		ChannelAdd = 4;
		NucContour = 1;		NucSegment = 2;		CellContour = 3;	CellSegment = 4;
	} // end of if, 2 segment option
	else {  // if only 1 segment type is chosen (e.g. color, B&W or Contour)
		if (LUToutline == "Yes") { ChannelAdd = 2;
			NucContour = 1;		CellContour = 2; 
		} else {
			if (LUTcolor == "Yes") { ChannelAdd = 2; }
			if (LUTbw == "Yes") { ChannelAdd = 2; }
			NucSegment = 1;		CellSegment = 2;
		}
	} // end of else, 2 segment option

	AmountIm = 0;	// sets to 0
	list = getFileList(filedir); // creates a list to go through all files in chosen folder
	setBatchMode(true);

	//	2  Open  ========================================================================================
	//	2.1  Open - Go through files in folder  ========================================================================================
	for (i=0; i<list.length; i++) { 	// go through all files in chosen folder
		showProgress(i+1, list.length);
		list = getFileList(filedir);
		print(i, " = ", list[i]);
		if(endsWith(list[i], "_montagePython_max.tif")){	 		// selects montage-max image
			ImageMax = list[i];
			AmountIm = AmountIm + 1;
			if (endsWith(ImageMax, delimiter2)) {	
				dotIndex = indexOf(ImageMax, delimiter2);
				LogFile = substring(ImageMax, 0, dotIndex);
			}
		}
		if(endsWith(list[i], "_montagePython_cell_mask.tif")){		 // selects montage-cell image
			ImageCell = list[i];
		}
		if(endsWith(list[i], "_montagePython_nucleus_mask.tif")){	 // selects montage-nucleus image
			ImageNuc = list[i];
		}
	}
	setBatchMode(false);

	//	2.2  Open - Open montage images  ========================================================================================
	print ("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print("* Amount of image sets to convert : ", AmountIm); 	// Amount of images (from an experiment to open)
	print("* File directory : ", filedir+ImageMax);			// File directory (chosen by user)
	print("* File images : ", LogFile);						// File name of directory
	open(filedir+ImageCell);
	rename(ImageCell);
	print("    - Cell segment image opened :        ", ImageCell);
	open(filedir+ImageNuc);
	rename(ImageNuc);
	print("    - Nucleus segment image opened : ", ImageNuc);
	open(filedir+ImageMax);
	rename(ImageMax);
	print("    - Max image opened :                       ", ImageMax);

	//	3  Adjust  ========================================================================================
	//	3.1  Adjust - Get dimensions  ========================================================================================
	selectWindow(ImageMax);
	getDimensions(widthMax, heightMax, channelsMax, slicesMax, framesMax);
	getPixelSize(unit, pixelWidth, pixelHeight);
	Column = widthMax / (width / Scale);
	Row = heightMax / (height / Scale);
	WithIm = widthMax / Column;
	HeightIm = heightMax / Row;
	print("* Amount of columns & rows : ", Column, " & " , Row, " (", WithIm, " * ", HeightIm," pixels)");

	//	3.1  Adjust - Create extra channels for segment images  ========================================================================================
	for (c = 1; c < ChannelAdd +1; c++) {
		Stack.setPosition(channelsMax, 1, 1);
		run("Add Slice", "add=channel");
		channelsM = channelsMax + 1;
	}
	print("* Amount of channels added : ", ChannelAdd);	
	if (LUToutline == "Yes") {
		if (LUTbwo == "No") {
			Stack.setPosition(channelsMax + CellContour, 1, 1);
			run("Add Slice", "add=channel");
			Stack.setPosition(channelsMax + CellContour + 1, 1, 1);
			print("       + Extra 'hidden' channel added, for cell detection ");
print("Cell segment image : ", channelsMax + CellContour + 1);
			selectWindow(ImageCell);
			run("Select All");
			run("Copy");
			selectWindow(ImageMax);
			Stack.setPosition(channelsMax + CellContour + 1, 1, 1);
			run("Paste");
		}
	}
	
	//	3.2  Adjust - Create  ========================================================================================
	print("* Type of LUT : ", LUT, "     [Outline-", LUToutline, ", Color-", LUTcolor, ", B&W-", LUTbw, ", Comibined-", LUTbwo, "]");

	//	3.2.1  Adjust - Create - Paste nuclear images into max projection  ========================================================================================
	if (LUTcolor == "Yes") {
		selectWindow(ImageNuc);
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + NucSegment, 1, 1);
		run("Paste");
		run("glasbey on dark");
		run("Enhance Contrast", "saturated=0.35");
		print("        - Nuclear segment : Glasbey on dark");
	}
	if (LUTbw == "Yes") {
		selectWindow(ImageNuc);
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + NucSegment, 1, 1);
		run("Paste");
		run("Grays");
		run("Enhance Contrast", "saturated=0.35");
		print("        - Nuclear segment : Black & white");
	}
	
	//	3.2.2  Adjust - Create - Paste cell images into max projection  ========================================================================================
	if (LUTcolor == "Yes") {
		selectWindow(ImageCell);
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + CellSegment, 1, 1);
		run("Paste");
		run("glasbey inverted");
		run("Enhance Contrast", "saturated=0.35");
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("00011");
		roiManager("Deselect");		roiManager("Show None");
	}
	if (LUTbw == "Yes") {
		selectWindow(ImageCell);
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + CellSegment, 1, 1);
		run("Paste");
		run("Grays");
		run("Enhance Contrast", "saturated=0.35");
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("00011");
		roiManager("Deselect");		roiManager("Show None");
	}

	//	3.2.3  Adjust - Create - Create nuclear contour  ========================================================================================
	selectWindow(ImageNuc);
	run("Select All");
	if (LUToutline == "Yes") {
		curROI = roiManager("count");
		if (curROI == 0) {
			roiManager("add");
		}
		ROImStart = roiManager("count") -1;
		getStatistics(areaSeg, meanSeg, minSeg, maxSeg, stdSeg, histogramSeg);
		for (m = 1; m < maxSeg+1; m++) {
			setThreshold(m, m);
			run("Create Selection");
			getStatistics(areaSeg2, meanSeg2, minSeg2, maxSeg2, stdSeg2, histogramSeg2);
			if (meanSeg2 == m) {
				roiManager("add");
				ROIm = roiManager("count") -1;
				roiManager("select", ROIm);
				run("Draw", "slice");
				roiManager("rename", "CellROI_"+m);
				resetThreshold();
			}
		}
		if (curROI == 0) {
			roiManager("select", 0);
			roiManager("Delete");
		}
		ROIselect = "CellROI_";		nR = roiManager("Count");	arrayROI_Cell = newArray();
		for (i=0; i<nR; i++) {		roiManager("Select", i);	rName = Roi.getName();
		if (indexOf(rName, ROIselect) >=0) { 	roiManager("Select", i);	arrayROI_Cell = Array.concat(arrayROI_Cell,i); 	}}
		if (LUToutline == "Yes") {
			Stack.setPosition(channelsMax + NucContour, 1, 1);
			roiManager("select", arrayROI_Cell);
			print("        - Nuclear outline : Magenta");
		}
		roiManager("select", arrayROI_Cell);
		roiManager("Combine");	roiManager("Add");	roiManager("delete");
		ROInuc = roiManager("count") -1;
		roiManager("select", ROInuc);
		roiManager("rename", "Nuc_Mask");
		roiManager("Set Color", "#FF00FF"); // magenta
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + NucContour, 1, 1);
		run("Paste");
		setMinAndMax(245, 250);
		run("Magenta");
	}

	//	3.2.4  Adjust - Create - Create cell contour  ========================================================================================
	selectWindow(ImageCell);
	run("Select All");
	if (LUToutline == "Yes") {
		ROImStart = roiManager("count") -1;
		getStatistics(areaSeg, meanSeg, minSeg, maxSeg, stdSeg, histogramSeg);
		for (m = 1; m < maxSeg+1; m++) {
			setThreshold(m, m);
			run("Create Selection");
			getStatistics(areaSeg2, meanSeg2, minSeg2, maxSeg2, stdSeg2, histogramSeg2);
			if (meanSeg2 == m) {
				roiManager("add");
				ROIm = roiManager("count") -1;
				roiManager("select", ROIm);
				run("Cyan");
				run("Draw", "slice");
				roiManager("rename", "CellROI_"+m);
				resetThreshold();
			}
		}
		ROIselect = "CellROI_";		nR = roiManager("Count");	arrayROI_Cell = newArray();
		for (i=0; i<nR; i++) {		roiManager("Select", i);	rName = Roi.getName();
		if (indexOf(rName, ROIselect) >=0) { 	roiManager("Select", i);	arrayROI_Cell = Array.concat(arrayROI_Cell,i); 	}}
		if (LUToutline == "Yes") {
			Stack.setPosition(channelsMax + CellContour, 1, 1);
			roiManager("select", arrayROI_Cell);
			print("        - Cellular outline : Cyan");
		}
		roiManager("select", arrayROI_Cell);
		roiManager("Combine");	roiManager("Add");	roiManager("delete");
		ROIcell = roiManager("count") -1;
		roiManager("select", ROIcell);
		roiManager("rename", "Cell_Mask");
		roiManager("Set Color", "#00FFFF"); // cyan
		run("Select All");
		run("Copy");
		selectWindow(ImageMax);
		Stack.setPosition(channelsMax + CellContour, 1, 1);
		run("Paste");
		setMinAndMax(245, 250);
		run("Cyan");
	}

	//	3.2.5  Adjust - Create - Show channels  ========================================================================================
	selectWindow(ImageMax);
	if (LUToutline == "Yes") {
		if (MaxChan == "Yes") {
			Stack.setActiveChannels("1110101");
			run("Brightness/Contrast...");
			for (i = 1; i < channelsMax +1; i++) {
				setSlice(i);
				run("Enhance Contrast", "saturated=0.35");
				getMinAndMax(min, max);
				max2 = max / 1.2;
				setMinAndMax(min, max2);
				print("        - Channel : ", i, " = min : ", min, " max : ", max2, " (1.5x of auto max)");
			}
		}
		if (MaxChan == "No") {
			Stack.setActiveChannels("0000101");
		}
		if (LUTbwo == "No") {
			if (MaxChan == "Yes") {
				Stack.setActiveChannels("111110");
			}
			if (MaxChan == "No") {
				Stack.setActiveChannels("000110");
			}
		}
	}
	if (LUTcolor == "Yes") {
		print("        - Cellular segment : Glasbey inverted");
		if (LUTbwo == "Yes") {
			Stack.setActiveChannels("0001111");
		}
		if (LUTbwo != "Yes") {
			Stack.setActiveChannels("00011");
		}
	}	
	if (LUTbw == "Yes") {
		print("        - Cellular segment : Black & white");
		if (LUTbwo == "Yes") {
			Stack.setActiveChannels("0001111");
		}
		if (LUTbwo != "Yes") {
			Stack.setActiveChannels("00011");
		}
	}
	roiManager("Deselect");		roiManager("Show None");

	//	3.2.6  Adjust - Create - Set position of montage image in screen  ========================================================================================
	getDimensions(WI, HE, CH, SL, FR);
	ratio = WI / HE;			ImageHeight = screenWidth * WI;			ImageScreenX = 0;				ImageScreenY = 20;
	Zoom = 1;
	if (Row > 5) {	Zoom = Zoom+1; ImageWidth = WI * 2; }	if (Row > 10) {	Zoom = Zoom+1; ImageWidth = WI * 1.5; }
	if (Row > 15) {	Zoom = Zoom+1; ImageWidth = WI * 1.5; }	if (Row > 20) {	Zoom = Zoom+1; ImageWidth = WI * 1.5; }
	setLocation(ImageScreenX, ImageScreenY, ImageWidth, screenHeight);
	for (z = 0; z < Zoom+1; z++) {
		run("In [+]");
	}
	run("Scale to Fit");
	print("       - Zoom fraction : ", d2s((getZoom()*100),2),"x");
	roiManager("Show None");
	roi = roiManager("count");
	if (roi != 0) {
		roiManager("delete");
	}
		
	//	4  Delete  ========================================================================================
	//	4.1  Delete - Dialog  ========================================================================================
	print ("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print ("     = = =    Delete Images / Cells    = = =");
	print ("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	getDimensions(widthM, heightM, channelsM, slicesM, framesM);
	Stack.setPosition(channelsM, 1, 1);
	
	Dialog.create("Delete Regions/Images ");
	Dialog.addMessage("> > >  Delete  < < <");
		items = newArray("Yes", "No");
	Dialog.addRadioButtonGroup("Delete images ? ", items, 1, 2, DeleteImages);
	Dialog.addNumber("Zoom image ? ", 1, 0, 4, " x  (for image selecting)");
	Dialog.addRadioButtonGroup("Delete cells/regions ? ", items, 1, 2, DeleteRegions);
	Dialog.addNumber("Zoom image ? ", 1, 0, 4, " x  (for cell/region selecting)");
	Dialog.addMessage("");
	Dialog.addMessage("");
	Dialog.show();
	// Delete
	DeleteImages = Dialog.getRadioButton();
	ZoomImage = Dialog.getNumber();
	DeleteRegions = Dialog.getRadioButton();
	ZoomCell = Dialog.getNumber();
	if (DeleteImages == "Yes") {
		ZoomCell = floor( ZoomCell / ZoomImage );
	}

	//	4.2  Delete - User interface  ========================================================================================
	//	4.2.1  Delete - User interface - Delete images  ========================================================================================
	if (DeleteImages == "Yes") {
		for (z = 1; z < ZoomImage; z++) {
			run("In [+]");
		}
		print("       - Zoom fraction : ", d2s((getZoom()*100),2),"x,  for clicking images to delete");
		setLocation(20, 20, screenWidth-100, screenHeight-40);
		selectWindow(ImageMax);
		ROIimagesStart = roiManager("count");
		setTool("point");
		run("Select None");
		run("Point Tool...", "type=Cross color=Pink size=[Extra Large] add label show");
		waitForUser(" Delete Images ", "    Select Images to delete    ");
		ROIimagesEnd = roiManager("count") -1;
		ROIcells = ROIimagesEnd - ROIimagesStart;
	}

	//	4.2  Delete - User interface  ========================================================================================
	//	4.2.2  Delete - User interface - Delete regions/cells  ========================================================================================
	if (DeleteRegions == "Yes") {
		for (z = 1; z < ZoomCell; z++) {
			run("In [+]");
		}
		print("       - Zoom fraction : ", d2s((getZoom()*100),2),"x,  for clicking cells/regions to delete");
		setLocation(20, 20, screenWidth-100, screenHeight-40);
		selectWindow(ImageMax);
		ROIcellsStart = roiManager("count");
		setTool("point");
		roiManager("Deselect");
		run("Point Tool...", "type=Dot color=Orange size=Large add label show");
		run("Select All");
		waitForUser(" Delete Regions / Cells ", "    Select Regions / Cells to delete    ");
		ROIcellsEnd = roiManager("count") -1;
		ROIcells = ROIcellsStart - ROIcellsStart;
	}

	//	4.3  Delete - Delete  ========================================================================================
	//	4.3.1  Delete - Delete - Delete images  ========================================================================================
	run("Scale to Fit");
	run("Out [-]");
	if (DeleteImages == "Yes") {
		for (c = ROIimagesStart; c < ROIimagesEnd +1; c++) {
			C = c - ROIimagesStart + 1;
			roiManager("select", c);
			Roi.getBounds(xCell,yCell,wCell,hCell);
			xCell = floor(xCell / WithIm) +1;
			yCell = floor(yCell / HeightIm) +1;
			Yname = newArray("000","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021");
			Xname = newArray("000","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021");
			CellString ="_" + Yname[yCell];
			CellString = CellString + "_" + Xname[xCell];
			roiManager("select", c);
			roiManager("rename", "DelImage_" + C + "_Image=" + Yname[yCell] + "_" + Xname[xCell]);
			print("        - Delete Image : ", C, "_Image=" + Yname[yCell] + "_" + Xname[xCell]);		
			
			// Open cell image to delete
			run("Image Sequence...", "open=&filedir file="+CellString+"_max_cell_mask.tif sort");
			CellSlice = getInfo("slice.label");
			if (endsWith(CellSlice, delimiter)) {	
				dotIndex = indexOf(CellSlice, delimiter);
				CellSlice = substring(CellSlice, 0, dotIndex);
			}
			rename(CellSlice);
			run("Select All");
			getStatistics(area, mean, min, max, std, histogram);
			if(min == max) {
				print("          ! Original image is already empty   ( Image=" + Yname[yCell] + "_" + Xname[xCell], " )");
				close(CellSlice);
			} else {
				run("Clear");
				saveAs("Tiff", filedir + CellSlice + ".tiff");
				rename(CellSlice);
				close(CellSlice);
				// Open nucleus image to delete
				run("Image Sequence...", "open=&filedir file="+CellString+"_max_nucleus_mask.tif sort");
				NucSlice = getInfo("slice.label");
				if (endsWith(NucSlice, delimiter)) {	
					dotIndex = indexOf(NucSlice, delimiter);
					NucSlice = substring(NucSlice, 0, dotIndex);
				}
				rename(NucSlice);
				run("Select All");
				run("Clear");
				saveAs("Tiff", filedir + NucSlice + ".tiff");
				rename(NucSlice);
				close(NucSlice);
			}
		}
	}

	//	4.3.2  Delete - Delete - Delete regions/cells  ========================================================================================
	if (DeleteRegions == "Yes") {
		for (c = ROIcellsStart; c < ROIcellsEnd +1; c++) {
			C = c - ROIcellsStart + 1;
			selectWindow(ImageMax);
			roiManager("select", c);
			getStatistics(areaAll, meanAll, minAll, maxAll, stdAll, histogramAll);
			Stack.setPosition(channelsM, 1, 1);
			Roi.getBounds(xCell,yCell,wCell,hCell);
			getRawStatistics(nPixelsCell, meanCell, minCell, maxCell, stdCell, histogramCell);
			xCell = floor(xCell / WithIm) +1;
			yCell = floor(yCell / HeightIm) +1;
			Yname = newArray("000","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021");
			Xname = newArray("000","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021");
			CellString ="_" + Yname[yCell];
			CellString = CellString + "_" + Xname[xCell];
			roiManager("select", c);
			roiManager("rename", "DelCell_" + C + "_Image=" + Yname[yCell] + "_" + Xname[xCell] + "_Cell=" + meanCell);
			print("        - Delete Cell : ", C, "_Image=" + Yname[yCell] + "_" + Xname[xCell] + "_Cell=" + meanCell);
			
			// Open cell image to delete
			run("Image Sequence...", "open=&filedir file="+CellString+"_max_cell_mask.tif sort");
			CellSlice = getInfo("slice.label");
			if (endsWith(CellSlice, delimiter)) {	// rename image if contains .tif
				dotIndex = indexOf(CellSlice, delimiter);
				CellSlice = substring(CellSlice, 0, dotIndex);
			}
			rename(CellSlice);						// end of rename
			getStatistics(area, mean, min, max, std, histogram);
			if(min == max) {
				print("          ! Original image doesn't contain cells", meanCell, "   ( Image=" + Yname[yCell] + "_" + Xname[xCell] + ", at Cell#=" + meanCell, " )");
				close(CellSlice);
			}
			if(min < max) {
			
				setThreshold(meanCell, meanCell);
				getStatistics(areaSel, meanSel, minSel, maxSel, stdSel, histogramSel);
				if (areaSel == areaAll) {
					print("          ! Original image doesn't contain cell", meanCell, "   ( Image=" + Yname[yCell] + "_" + Xname[xCell] + ", at Cell#=" + meanCell, " )");
				}
				if (areaSel < areaAll) {
					run("Create Selection");
					roiManager("add");
					roiManager("Show None");
					resetThreshold();
					ROIcellSelect = roiManager("count") -1;
					roiManager("select", ROIcellSelect);
					roiManager("rename", "DelCell_Cell_" + meanCell);
					roiManager("select", ROIcellSelect);
					run("Clear");
					saveAs("Tiff", filedir + CellSlice + ".tiff");
					rename(CellSlice);
					close(CellSlice);
					// Open nucleus image to delete
					run("Image Sequence...", "open=&filedir file="+CellString+"_max_nucleus_mask.tif sort");
					NucSlice = getInfo("slice.label");
					if (endsWith(NucSlice, delimiter)) {	
						dotIndex = indexOf(NucSlice, delimiter);
						NucSlice = substring(NucSlice, 0, dotIndex);
					}
					rename(NucSlice);
					roiManager("select", ROIcellSelect);
					run("Clear");
					saveAs("Tiff", filedir + NucSlice + ".tiff");
					rename(NucSlice);
					close(NucSlice);
				}
			}
			if(min >= max) {
				close(CellSlice + ".tif");
			}
		}
	}

	//	5  Finish  ========================================================================================
	//	5.1  Finish - Save  ========================================================================================
	roiManager("Deselect");
	roiManager("List");
	Table.rename("Overlay Elements", "Results_"+LogFile);
	print ("* ROI/Log Text File Saved : " + filedir + "/" + LogFile  + "_montageROIdelete_["+MACRO+"_"+DateString+"].txt");
	saveAs("Results", filedir + "/" + LogFile  + "_montageROIdelete_["+MACRO+"_"+DateString+"].txt");
	Table.rename(LogFile  + "_montageROIdelete_["+MACRO+"_"+DateString+"].txt", "Results_"+LogFile);
	
	//	5.2  Finish - Close  ========================================================================================
	Dialog.create("Close Images");
		items = newArray("Yes", "No");
	Dialog.addRadioButtonGroup("Ready to close all images / data ? ", items, 1, 2, "Yes");
	Dialog.show();
	Close = Dialog.getRadioButton();
	EndTime = d2s((((getTime()-StartTime)/1000)/60),2);	
	print("* Time spent on deleting images/regions : ", EndTime, " min");
	if (Close == "Yes") {
		run("Close All");																			// close all open images
		if (isOpen("Results_"+LogFile)) {	selectWindow("Results_"+LogFile);	run("Close");	}	// close results
		if (isOpen("ROI Manager")) {		selectWindow("ROI Manager");		run("Close");	}	// close ROI manager
		if (isOpen("Log")) {				selectWindow("Log");				
		saveAs("Text", filedir + "/" + LogFile  + "_montageLog_["+MACRO+"_"+DateString+"].txt");	run("Close");	}	// close Log file
	}
	
	//	5.3  Finish - End of macro  ========================================================================================
}