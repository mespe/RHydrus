
hydrusSA = 
function(args)
{    
   .Call("R_hydrus_standalone", as.character(c("hydrus", args)))
}
