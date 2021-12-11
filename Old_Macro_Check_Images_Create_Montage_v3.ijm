/*   Macro to check images for bubbles, intensity, etc

Version:  See below
Description: 
   * This macro can create a montage image (or animation) of all images in one folder with the option to open from different extensions.
   * You have the option to save and/or copy the image at the end
   * 
To do:
   * Create option to delete regions (in mask images) or a list of which images to delete (in logfile)

Contact:   Linda Joosen
*/

//  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  
// Variables 
	Ch1 = "Cy3";
	Ch2 = "Cy5";
	Ch3 = "Dapi";
	DefContrast = "No"; 	// default setting for dialog - Change contrast
	DefCopy = "Yes";		// default setting for dialog - Copy montage image
	DefSave = "No";		// default setting for dialog - Save montage image
	
//  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  ///  //  

// Macro credentials
	title		=	"Macro_Check_Images_Create_Montage";
	version		=	"v3";
	MACRO 		= 	title+version;
	date		=	"24 Oct 2019";
	Contact		= 	"l.joosen@nki.nl"
	Internet	= 	"https://spark.adobe.com/page/xKqkGsRhW63IS/";
// Variables (changing will affect the script)
	ImageType1 = "_max.tif";					channel = 3;	c1 = Ch1;	c2 = Ch2;		c3 = Ch3;
	ImageType2 = "_max_nucleus_mask.tif";		channel = 1;	c1 = "_nucleus_mask";
	ImageType3 = "_max_cell_mask.tif";			channel = 1;	c1 = "_cell_mask";
	ImageType4 = "_mask_cell+nuc+spots.tif";	channel = 6;	c1 = Ch1;	c2 = Ch1+"_";	c3 = Ch1+"_";	c4 = Ch1+"_";	c5 = Ch1+"_";	c6 = Ch1+"_";
	ImageType5 = "_loc_results_cy3.tif";		channel = 6;	c1 = Ch1;	c2 = Ch1+"_";	c3 = Ch1+"_";	c4 = Ch1+"_";	c5 = Ch1+"_";	c6 = Ch1+"_";
	ImageType6 = "_loc_results_cy5.tif";		channel = 6;	c1 = Ch2;	c2 = Ch2+"_";	c3 = Ch2+"_";	c4 = Ch2+"_";	c5 = Ch2+"_";	c6 = Ch2+"_";
	ImageType7 = "+TS.tif";						channel = 4;	c1 = Ch1;	c2 = Ch2;		c3 = Ch1+"_TS";	c4 = Ch2+"_TS";
	Scale = 0.25;

macro title {
// Start dialog
	Dialog.create("Open Images");
	Dialog.addMessage("> > >  Open  < < <");
	Dialog.addChoice("Work on Image : ", newArray("Current file", "Open sequence"), "Open sequence");
	Dialog.addChoice("Images contain : ", newArray("N/A", ImageType1, ImageType2, ImageType3, ImageType4, ImageType5, ImageType6, ImageType7), ImageType1);
	//Dialog.addMessage("");
	//Dialog.addMessage("> > >  Adjustment  < < <");
	//Dialog.addChoice("Create hyperstack : ", newArray("Yes", "No"), "Yes");
	//Dialog.addToSameRow();
	//Dialog.addMessage(" (only if current selected)");
	Dialog.addMessage("");
	Dialog.addMessage("> > >  Create  < < < ");
	Dialog.addChoice("Create montage : ", newArray("Not", "Manual", "Automatic"), "Automatic");
	Dialog.addNumber("Animation speed : ", 1.500, 1, 7.5," (if montage isn't performed) ");
	Dialog.show();
	WorkOn = Dialog.getChoice();
	Contain = Dialog.getChoice();
	//CreateHyperstack = Dialog.getChoice();
	CreateMontage = Dialog.getChoice();
	AnimationSpeed =  Dialog.getNumber();

//===> 0.4    <========= Print Date =======================================================
	print ("~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~~=>+<=~"); 
	//	1.3.1  Date & time:
	print("+~=-=~+");
	selectWindow("Log");
	print("\\Clear");
	print("Macro : ", title + ": " + version); 

	//	1.3.1  Date & time:
	MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
		TimeString ="Date: "+DayNames[dayOfWeek]+" ";
		if (dayOfMonth<10) {TimeString = TimeString+"0";}
		TimeString = TimeString+dayOfMonth+" "+MonthNames[month]+", "+year+"\nTime: ";
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
	print("Dimensions desktop screen : ", screenWidth, " * ", screenHeight);
	print("Contact : ",Contact );
	print("Info : see line below    V ");
	print(Internet );
	//print("Contact : ",Contact ,",   for more information  (check link in next line)");
	//print(Internet);
	print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");

//===> 1      <========= Open Files =======================================================
	print("1 # Image information :");
	print ("       * Image used :            " + WorkOn); 	
	if (WorkOn == "Current file" ) {
		filedir = getDirectory("image");
		FileName = getTitle();
		getDimensions(Width, Height, Channels, Slices, Frames);
	}
	if (WorkOn != "Current file" ) {
		run("Close All");
		filedir = getDirectory("Choose a Directory");
		list = getFileList(filedir);
		im = list.length / channel;
		print ("       * Images consists :      " + Contain);
		print ("       * Amount of images :   " + list.length);
		if(Contain != "N/A"){
			run("Image Sequence...", "open=&filedir file="+Contain+" sort");
		}
		if(Contain == "N/A"){
			run("Image Sequence...", "open=&filedir sort");
		}
		FileName = getTitle();
		getDimensions(Width, Height, Channels, Slices, Frames);
		rename(FileName);
		filein=filedir+FileName;
	}
	print ("       * File name : ", FileName);
	print ("       * Directory : ", filedir);

	// Define amount of channels
	getDimensions(width, height, channels, slices, frames);
	if ( Contain == ImageType1 ) {	channel = 3;	c1 = Ch1;	c2 = Ch2;		c3 = Ch3;
		Ch1Color = "Red";	Ch2Color = "Green";	Ch3Color = "Blue";	Ch4Color = "Grays";	Ch5Color = "Cyan";	Ch6Color = "Magenta";	}
	if ( Contain == ImageType2 ) {	channel = 1;	c1 = "_nucleus_mask";
		Ch1Color = "Grays";	}
	if ( Contain == ImageType3 ) {	channel = 1;	c1 = "_cell_mask";
		Ch1Color = "Grays";	}
	if ( Contain == ImageType4 ) {	channel = 6;	c1 = Ch1;	c2 = Ch1+"_";	c3 = Ch1+"_";	c4 = Ch1+"_";	c5 = Ch1+"_";	c6 = Ch1+"_";
		Ch1Color = "Red";	Ch2Color = "Green";	Ch3Color = "Blue";	Ch4Color = "Grays";	Ch5Color = "Cyan";	Ch6Color = "Magenta";	}
	if ( Contain == ImageType5 ) {	channel = 6;	c1 = Ch1;	c2 = Ch1+"_";	c3 = Ch1+"_";	c4 = Ch1+"_";	c5 = Ch1+"_";	c6 = Ch1+"_";
		Ch1Color = "Red";	Ch2Color = "Green";	Ch3Color = "Blue";	Ch4Color = "Grays";	Ch5Color = "Cyan";	Ch6Color = "Magenta";	}
	if ( Contain == ImageType6 ) {	channel = 6;	c1 = Ch2;	c2 = Ch2+"_";	c3 = Ch2+"_";	c4 = Ch2+"_";	c5 = Ch2+"_";	c6 = Ch2+"_";
		Ch1Color = "Red";	Ch2Color = "Green";	Ch3Color = "Blue";	Ch4Color = "Grays";	Ch5Color = "Cyan";	Ch6Color = "Magenta";	}
	if ( Contain == ImageType7 ) {	channel = 4;	c1 = Ch1;	c2 = Ch2;		c3 = Ch1+"_TS";	c4 = Ch2+"_TS";
		Ch1Color = "Red";	Ch2Color = "Green";	Ch3Color = "Blue";	Ch4Color = "Grays";	Ch5Color = "Cyan";	Ch6Color = "Magenta";	}

// Set contrast & color of each channel
	run("Out [-]");
	getDimensions(width, height, channels, slices, frames);
	AllSlices = channels * slices * frames;
	if (AllSlices == slices) { 
		CreateHyperstack = "Yes";
	}
	
	if (CreateHyperstack == "Yes") {
			Frames = slices / channel;
			print ("       * Dimensions of image :            " + Frames, slices, channel); 	
	} else {
		Frames = frames;
	}
	print ("       * Filename : ", FileName);
	print ("       * Amount of channels : ", channel, ", Frames : ", Frames);

	if (CreateHyperstack == "Yes") {
		getDimensions(width, height, channels, slices, frames);	
		run("Stack to Hyperstack...", "order=xyczt(default) channels="+channel+" slices=1 frames="+Frames+" display=Color");
	}
	if (frames > 2) {
		FramesHalf = frames / 2;
	}
	for (c = 1; c < channel+1; c++) {
		Stack.setPosition(c, 1, 1);
		run("Brightness/Contrast...");	run("Enhance Contrast", "saturated=0.35");	
		getMinAndMax(min, max);			
		if (frames > 1) { 
			Stack.setPosition(c, 1, frames);
			run("Brightness/Contrast...");	run("Enhance Contrast", "saturated=0.35");	
			getMinAndMax(min2, max2);
			if (min2 < min) { min = min2; } else { min = min; } 	if (max2 > max) { 	max = max2; } else { max = max; }
		}
		if (frames > 2) { 
			Stack.setPosition(c, 1, FramesHalf);
			run("Brightness/Contrast...");	run("Enhance Contrast", "saturated=0.35");	
			getMinAndMax(min2, max2);			
			if (min2 < min) { min = min2; } else { min = min; } 	if (max2 > max) { 	max = max2; } else { max = max; }
		}
		print ("       * Channel ", c, " : min-max : ", min, "-" , max);
		if ( Contain != ImageType4 ) {
			if (c == 1) { min1 = minA = min; max1 = maxA = max; run(Ch1Color); }	if (c == 2) { min2 = minA = min; max2 = maxA = max; run(Ch2Color); }	if (c == 3) { min3 = minA = min; max3 = maxA = max; run(Ch3Color); }
			if (c == 4) { min4 = minA = min; max4 = maxA = max; run(Ch4Color); }	if (c == 5) { min5 = minA = min; max5 = maxA = max; run(Ch5Color); }	if (c == 6) { min6 = minA = min; max6 = maxA = max; run(Ch6Color);	}
		}		
		if ( Contain == ImageType4 ) { 
			if (c == 1) { minA = min; maxA = max * 1.6; run(Ch1Color); }	if (c == 2) { minA = min; max = maxA * 2.5; run(Ch2Color); }	if (c == 3) { minA = min; maxA = max * 3; run(Ch3Color); }
			if (c == 4) { minA = 0; maxA = 60; run(Ch4Color); }				if (c == 5) { minA = 0; maxA = 1; run(Ch5Color); }				if (c == 6) { minA = 0; maxA = 1; run(Ch6Color);	}
		}
		setMinAndMax(minA, maxA);
		run("Channels Tool...");
	}
	if ( Contain == ImageType4 ) { 
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("100110");
		run("In [+]");	run("In [+]");	run("In [+]");	run("In [+]");
	}
	Row = d2s(Frames / 3, 0);
	Column = Frames / Row;	Column = floor(Column+1);
	Row = Frames / Column;	Row = floor(Row+1);
	getDimensions(width, height, channels, slices, frames);
	if (CreateMontage != "Not") {
		run("Rotate 90 Degrees Right");
		run("Flip Horizontally", "stack");
		print ("       * Note: Image is rotated 90 degrees right & horizontally flipped, for creating fitting montage image");


	print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print("2 # Montage information :");
// create montage image by automatic calculation			
		if (CreateMontage == "Automatic") {
			if (Contain == ImageType1) {	SliceMin = 3; 	}	if (Contain == ImageType2) {	SliceMin = 5; 	}	if (Contain == ImageType3) {	SliceMin = 5; 	}	
			if (Contain == ImageType4) {	SliceMin = 4;	}	if (Contain == ImageType5) {	SliceMin = 5; 	}	if (Contain == ImageType6) {	SliceMin = 5; 	}	
			if (Contain == ImageType7) {	SliceMin = 2; 	}	
			for (fr = 1; fr < frames+1; fr++) {
				Stack.setPosition(1, 1, fr);
				SliceInfo = getInfo("slice.label");
				SliceSplit = split(SliceInfo, "_");
				Co = SliceSplit.length - SliceMin+1;
				Ro = SliceSplit.length - SliceMin;
				for(f=1; f<SliceSplit.length+1; f++) {
					if (f == Co) {	Column = SliceSplit[f]; 	Column = parseInt(Column);	}
					if (f == Ro) {	Row = SliceSplit[f]; 		Row = parseInt(Row);	}
				}
				if (fr == 1) {	Row1 = Row; }
				if (fr != 1) {	if (Row > Row1) { Row1 = Row; } }
			}
				Column = Column+1;		Row = Row1 +1;
		}
// create montage image by defining column
		if (CreateMontage == "Manual") {
			Dialog.create("Montage");
			Dialog.addMessage("");
			Dialog.addMessage("Total frames : "+ Frames);
			Dialog.addMessage("");
			Dialog.addNumber("Columns : ", Column);
			// Dialog.addNumber("Rows : ", Row);
			Dialog.show();
			Column = Dialog.getNumber();
			// Row = Dialog.getNumber();
			Row = frames / Column; 		Row = floor(Row);
		}
		print ("       * Total images : ", Frames, " = Row(s): ", Row, ", Column(s): ", Column);
// Text settings
		setColor(255, 255, 255);
		FontSize = height /150;
		setFont("Arial", FontSize);
	
// Create Montage image from sequence images
			run("Make Montage...", "columns="+Column+" rows="+Row+" scale="+Scale+" first=1 last="+Frames+" increment=1 border=0 font="+FontSize+" label");
// Set contrast for segmentation images
		if ( Contain == ImageType2) {
			setMinAndMax(0, 0.5);
		}
		if ( Contain == ImageType3) {
			setMinAndMax(0, 0.5);
		}
// Set size of montage image to fit screen
		getDimensions(widthM, heightM, channelsM, slicesM, framesM);
		ratio = width / height;
		ImageHeight = screenHeight / 0.9;			ImageWidth = ImageHeight * ratio;
		ImageScreenX = 0;							ImageScreenY = 20;
		setLocation(ImageScreenX, ImageScreenY, ImageWidth, ImageHeight);

		Dialog.create("Delete & Save");
		Dialog.addMessage("> > >  Delete  < < <");
		Dialog.addChoice("Delete regions ? ", newArray("Yes", "No"), "No");
		Dialog.addChoice("Delete Images ? ", newArray("Yes", "No"), "No");
		if (Contain == ImageType1) {	item = 3;	item1 = c2;		items = newArray(c1, c2, c3);	c4 = "-";	c5 = "-";	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 1, item, item1);	}
		if (Contain == ImageType2) {	item = 1;	item1 = c1;		items = newArray(c1);						c2 = "-";	c3 = "-";	c4 = "-";	c5 = "-";	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 1, item, item1);	}
		if (Contain == ImageType3) {	item = 1;	item1 = c1;		items = newArray(c1);						c2 = "-";	c3 = "-";	c4 = "-";	c5 = "-";	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 1, item, item1);	}
		if (Contain == ImageType4) {	item = 6;	item1 = c1;		items = newArray(c1, c2, c3, c4, c5, c5);	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 2, item/2, item1);	}
		if (Contain == ImageType5) {	item = 6;	item1 = c1;		items = newArray(c1, c2, c3, c4, c5, c5);	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 2, item/2, item1);	}
		if (Contain == ImageType6) {	item = 6;	item1 = c1;		items = newArray(c1, c2, c3, c4, c5, c5);	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 2, item/2, item1);	}
		if (Contain == ImageType7) {	item = 4;	item1 = c2;		items = newArray(c1, c2, c3, c4);			c5 = "-";	c6 = "-";
			Dialog.addRadioButtonGroup("  Select from channel :                ", items, 1, item, item1);	}
		Dialog.addMessage("");
		Dialog.addMessage("> > >  Copy/Save  < < <");
		//Dialog.addNumber("Scale montage image ? ", 1, 0, 7.5," (1 = not scaling, 2 = 2x , etc) ");
		Dialog.addChoice("Change contrast manually ? ", newArray("Yes", "No"), DefContrast);
		Dialog.addChoice("Contrast type ? ", newArray("Cell value", "All same"), "All same");
		Dialog.addToSameRow();
		Dialog.addMessage("= Only for mask images ");
		Dialog.addChoice("Copy montage image ? ", newArray("Yes", "No"), DefCopy);
		Dialog.addChoice("Save montage image ? ", newArray("Yes", "No"), DefSave);
		Dialog.show();
		DeleteRegions = Dialog.getChoice();
		DeleteImages = Dialog.getChoice();
		DeleteChan = Dialog.getRadioButton;
		//ScaleImage = Dialog.getNumber();
		ContrastImage = Dialog.getChoice();
		ContrastCells = Dialog.getChoice();
		CopyImage = Dialog.getChoice();
		SaveImage = Dialog.getChoice();
		if (DeleteChan == c1) { chan = 1; }	if (DeleteChan == c2) { chan = 2; }	if (DeleteChan == c3) { chan = 3; }
		if (DeleteChan == c4) { chan = 4; }	if (DeleteChan == c5) { chan = 5; }	if (DeleteChan == c6) { chan = 6; }

		Stack.setPosition(chan, 1, 1);
		//if (ScaleImage  > 1) {
		//	widthMscale = widthM / ScaleImage; 	heightMscale = heightM / ScaleImage;
		//	run("Scale...", "x="+widthMscale+" y="+heightMscale+" z="+slicesM); 
		//}

	print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print("3 # Copy & save information :");

		// creates user define contrast option
		if (ContrastImage  == "Yes") {
			getDimensions(wid, hei, chan, sli, fra);
			print ("       * Contrast images are manually adjusted, for ", chan, " channels ");	
			run("Brightness/Contrast...");	COL = "No";
			if ( Contain != ImageType2 ) {	COL = "Yes"; }
			if ( Contain != ImageType3 ) {	COL = "Yes"; }
			if ( COL == "Yes" ) {	Stack.setDisplayMode("color"); }
			for (c = 1; c < chan+1; c++) {
				Stack.setPosition(c, 1, 1);
				waitForUser("Change contast channel : " + c);
				getMinAndMax(min, max);
				print ("       * Contrast (min-max) channel ", c, " = ", min, "-", max);	
			}
			if ( COL == "Yes" ) {	Stack.setDisplayMode("composite"); }
		}
		// set contrast of segmented images to single colors
		if (ContrastCells  != "All same") {		CON = "No";
			if ( Contain == ImageType2 ) {	CON = "Yes"; }
			if ( Contain == ImageType3 ) {	CON = "Yes"; }
			if ( CON == "Yes" ) {	run("Enhance Contrast", "saturated=0.35"); }
		}
		// save montage image
		if (SaveImage  == "Yes") {
			saveAs("Tiff", filedir+FileName+"_Montage");
			print("Montage tiff image is saved as : ", filedir+FileName);
		}
		// copy montage image
		if (CopyImage  == "Yes") {
			run("Copy to System");
			waitForUser("Your " + Contain + " montage image is copied ");
		}

		Del = "No";
		if (DeleteRegions == "Yes") { Del = "Yes"; }	if (DeleteImages == "Yes") { Del = "Yes"; }
		if (Del == "Yes") {

	print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
	print("4 # Delete information :");
	
			if (DeleteRegions == "Yes") {
				ROI = roiManager("count");
				if (ROI > 0) {
					roiManager("deselect");	roiManager("delete");
				}
				selectWindow("Montage");
				getDimensions(widthM, heightM, channelsM, slicesM, framesM);
				widthMs = widthM / Scale;	heightMs = heightM / Scale;
				setTool("wand");
				ROIstart = roiManager("count");
				roiManager("Show All with labels");
				waitForUser("DELETE REGIONS/CELLS ", "Click on regions to be deleted, confirm each with Ctrl/Cmnd + T");
				ROIend = roiManager("count");
				ROIsel = ROIend - ROIstart;
				print ("       * ROIs to delete : ", ROIstart, " - " , ROIend, " (", ROIsel, "x)");
				for (r = 0; r < ROIsel; r++) {
					roiManager("select", r);
					getBoundingRect(xa, ya, wa, ha);
					xaS = xa / Scale;	yaS = ya / Scale;	waS = wa / Scale;	haS = ha / Scale;
					ImageH = yaS / height; 	ImageH = floor(ImageH);
					ImageW = xaS / width; 	ImageW = floor(ImageW)+1;
					ImH = ImageH;	ImW = ImageW - 1;
					xaS = (xaS - (width * ImW));	yaS = (yaS - (height * ImH));
					Imageslice = ImageH * Column;
					ImageSlice = Imageslice + ImageW;
					selectWindow(FileName);
					setSlice(ImageSlice);
					makeRectangle(xaS, yaS, waS, haS);
					roiManager("add");
					roiSel = roiManager("count") - 1;
					roiManager("select", roiSel);
					rPlus = r + 1;
					roiManager("rename", "R"+rPlus+"_S"+ImageSlice+"_(x"+xaS+"-y"+yaS+"-w"+waS+"-h"+haS+")");
				}
			} // end of delete regions
			if (DeleteImages == "Yes") {
				selectWindow("Montage");
				getDimensions(widthM, heightM, channelsM, slicesM, framesM);
				widthMs = widthM / Column; heightMs = heightM / Row;
				for (w = 1; w < Column+1; w++) {
					W = w * widthMs;
					setForegroundColor(250, 250, 250);
					drawLine(W, 0, W, heightM);
				}
				for (h = 1; h < Row+1; h++) {
					H = h * heightMs;
					setForegroundColor(250, 250, 250);
					drawLine(0, H, widthM, H);
				}
				ROIstart = roiManager("count");
				setTool("point");
				run("Point Tool...", "type=Circle color=Yellow size=[Extra Large] add label show");
				roiManager("Show All with labels");
				waitForUser("Images to delete : ", "Click in images to be deleted. Press OK when ready");
				ROIend = roiManager("count");
				ROIsel = ROIend - ROIstart;
				print ("       * Images to delete : ", ROIstart, " - " , ROIend, " (", ROIsel, "x)");
				arrayROI = newArray("TO DELETE : ");
				for (r = 0; r < ROIsel; r++) {
					roiManager("select", ROIstart+r);
					getBoundingRect(xa, ya, wa, ha);
					xaS = xa / Scale;	yaS = ya / Scale;
					ImageH = yaS / height; 	ImageH = floor(ImageH);
					ImageW = xaS / width; 	ImageW = floor(ImageW)+1;
					ImH = ImageH;	ImW = ImageW - 1;
					xaS = (xaS - (width * ImW));	yaS = (yaS - (height * ImH));
					Imageslice = ImageH * Column;
					ImageSlice = Imageslice + ImageW;
					selectWindow(FileName);
					setSlice(ImageSlice);
					SliceInfo = getInfo("slice.label");			
					arrayROI = Array.concat(arrayROI,SliceInfo);
				}
				roiManager("select", arrayROI);
				Array.print(arrayROI);
			} // end of delete images
		} // end of Delete = all
	} // end of montage

	if (CreateMontage == "Not") {
		run("Animation Options...", "speed="+AnimationSpeed+" start");
		}

	print ("       * Montage image : ", Frames, " (columns = ", Column, " & rows = ", Row);

	Dialog.create("Close");
	Dialog.addMessage("> > >  Close  < < <");
	Dialog.addMessage("");
	Dialog.addMessage("File name : " + FileName);
	Dialog.addMessage("");
		items = newArray("Yes", "No");
	Dialog.addRadioButtonGroup("   Ready to close all images ? ", items, 1, 2, "Yes");
	Dialog.addMessage("");
	Dialog.addMessage("    - Total frames : " + Frames);
	Dialog.addMessage("    - Column(s) : " + Column);
	Dialog.addMessage("    - Row(s) : " + Row);
	Dialog.show();
	CloseAll = Dialog.getRadioButton;

	if (CloseAll == "Yes") {
	 run("Close All");
	}
	
} // end of macro all

	