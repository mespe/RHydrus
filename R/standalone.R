
hydrusSA = 
function(args = character())
{    
   .Call("R_hydrus_standalone", as.character(c("hydrus", args)))
}
