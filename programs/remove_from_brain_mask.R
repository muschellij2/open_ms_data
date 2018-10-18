x = list.files(pattern = ".nii.gz", 
	recursive=TRUE)
xx = x[ grepl("prenorm", x)]
removers = grepl(
	paste0("(B|b)rain|WM|CSF|TISSUE|STRUCT", 
		"|FAST|resampled|regtoMNI|regtoEve"), 
	basename(xx))
brain = xx[removers]
unique(basename(brain))
unique(basename(xx[!removers]))


file.remove(brain)


x = list.files(
	pattern = ".nii.gz", 
	recursive=TRUE)
xx = x[ grepl("malf/", x)]
file.remove(xx)


x = list.files(
	pattern = ".mat$", 
	recursive=TRUE)
xx = x[ grepl("malf/", x)]
file.remove(xx)
