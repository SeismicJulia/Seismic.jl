function SeisCopy(in,out)

	run(`cp ${join([in ".seisd"])} ${join([out ".seisd"])}`);run(`cp ${join([in ".seish"])} ${join([out ".seish"])}`);
	
end
