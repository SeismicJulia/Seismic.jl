function SeisRemove(in::ASCIIString)
	rm(join([in ".seisd"]));rm(join([in ".seish"]));
end