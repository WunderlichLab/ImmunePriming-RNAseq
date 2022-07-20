macro "Save Images [1]"{
	image_dump = getDirectory("Image Storage")
	i = 1
	j = 1
	waitForUser;
	while (i < 7){
		run("Duplicate...", " ");
		title = getTitle();
		saveAs("Tiff", image_dump + i + "_" + j + ".tif");
		close();
		j = j + 1;
		if (j == 4) {
			i = i + 1;
			j = 1;
		}
		waitForUser;
	}
}

macro "Count 15-pix [5]"{
	input_path = getDirectory("input files");
	output_path = input_path
	fileList = getFileList(input_path)
	
	for (i=0; i<fileList.length; i++){
	
		//Clean-up to prepare for next image
		roiManager("reset");
		run("Close All");
		run("Clear Results");
	
		//Open image
		open(input_path + fileList[i]);
		print(input_path + fileList[i]); //displays file that is processed
		title = getTitle();
	
		//Coding Block
		run("16-bit");
		run("Threshold...");
		waitForUser;
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");
		run("Watershed");
		run("Analyze Particles...", "size=15-Infinity circularity=0.2-1.00 exclude clear add");
	
		saveAs("results", output_path+title+"_results.txt");
	
	}
}

macro "Count 30-pix [6]"{
	input_path = getDirectory("input files");
	output_path = input_path
	fileList = getFileList(input_path)
	
	for (i=0; i<fileList.length; i++){
	
		//Clean-up to prepare for next image
		roiManager("reset");
		run("Close All");
		run("Clear Results");
	
		//Open image
		open(input_path + fileList[i]);
		print(input_path + fileList[i]); //displays file that is processed
		title = getTitle();
	
		//Coding Block
		run("16-bit");
		run("Threshold...");
		waitForUser;
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");
		run("Watershed");
		run("Analyze Particles...", "size=30-Infinity circularity=0.2-1.00 exclude clear add");
	
		saveAs("results", output_path+title+"_results.txt");
	
	}
}

macro "Count 50-pix [7]"{
	input_path = getDirectory("input files");
	output_path = input_path
	fileList = getFileList(input_path)
	
	for (i=0; i<fileList.length; i++){
	
		//Clean-up to prepare for next image
		roiManager("reset");
		run("Close All");
		run("Clear Results");
	
		//Open image
		open(input_path + fileList[i]);
		print(input_path + fileList[i]); //displays file that is processed
		title = getTitle();
	
		//Coding Block
		run("16-bit");
		run("Threshold...");
		waitForUser;
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");
		run("Watershed");
		run("Analyze Particles...", "size=50-Infinity circularity=0.2-1.00 exclude clear add");
	
		saveAs("results", output_path+title+"_results.txt");
	
	}
}

macro "Count 80-pix [8]"{
	input_path = getDirectory("input files");
	output_path = input_path
	fileList = getFileList(input_path)
	
	for (i=0; i<fileList.length; i++){
	
		//Clean-up to prepare for next image
		roiManager("reset");
		run("Close All");
		run("Clear Results");
	
		//Open image
		open(input_path + fileList[i]);
		print(input_path + fileList[i]); //displays file that is processed
		title = getTitle();
	
		//Coding Block
		run("16-bit");
		run("Threshold...");
		waitForUser;
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");
		run("Watershed");
		run("Analyze Particles...", "size=80-Infinity circularity=0.2-1.00 exclude clear add");
	
		saveAs("results", output_path+title+"_results.txt");
	
	}
}

macro "Count 120-pix [9]"{
	input_path = getDirectory("input files");
	output_path = input_path
	fileList = getFileList(input_path)
	
	for (i=0; i<fileList.length; i++){
	
		//Clean-up to prepare for next image
		roiManager("reset");
		run("Close All");
		run("Clear Results");
	
		//Open image
		open(input_path + fileList[i]);
		print(input_path + fileList[i]); //displays file that is processed
		title = getTitle();
	
		//Coding Block
		run("16-bit");
		run("Threshold...");
		waitForUser;
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Watershed");
		run("Watershed");
		run("Analyze Particles...", "size=120-Infinity circularity=0.2-1.00 exclude clear add");
	
		saveAs("results", output_path+title+"_results.txt");
	
	}
}