args <- commandArgs(T)
for (i in 1:length(args)) {
	if (i==1) { instruction <- paste(args[i]) }
	if (i>1) instruction <- paste(instruction, args[i])
}
print( instruction )
system( instruction )
