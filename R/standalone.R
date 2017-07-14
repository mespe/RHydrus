
hydrusSA = 
function(args = character())
{    
   .Call("R_hydrus_standalone", as.character(c("hydrus", args)), PACKAGE = "RHydrus")
}



hydrusF = 
function(options = NA, selector = NA, profile = NA, atmosph = NA,
         profile_data = NULL,
         indir = ".", outdir = ".", showOutput = FALSE, args = character(),
         callbackFun = NULL)
{
   files = c(options = options, selector = selector, profile = profile, atmosph = atmosph)
   if(!missing(indir) && any(is.na(files))) {
     files = getInputFilesFromDir(indir, files)
   } else {
     files = path.expand(files)
   }
   

   if(any(w <- is.na(files))) {
      warning("missing input files: ", paste(files[w], collapse = ", "))
      files[w] = ""
   }
      
   if(!file.exists(outdir) && !dir.create(outdir))
      stop("can't create output directory ", outdir)
       
   if(length(profile_data))
      .Call("R_hydrus_direct", files["options"], files["selector"], files["atmosph"], files["profile"], outdir, as.logical(showOutput), profile_data, PACKAGE = "RHydrus")       
   else
      .Call("R_hydrus_inputs", files["options"], files["selector"], files["atmosph"], files["profile"], outdir, as.logical(showOutput), callbackFun, PACKAGE = "RHydrus")
}

getInputFilesFromDir =
function(dir, files, ...)
{
   w = is.na(files)
   if(any(w)) {
       dirFiles = list.files(dir, full = TRUE)
       files[w] = sapply(names(files[w]), findFile, dirFiles, ...)
   }

   structure(path.expand(files), names = names(files))
}

findFile =
function(pattern, files, ...)
{
   i = grep(pattern, basename(files), ignore.case = TRUE, ...)
   if(length(i))
      files[i][1]
   else
      NA
}
