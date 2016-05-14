if [ $NERSC_HOST == "edison" ]; then
	nsockets=2
	ncore_per_socket=12
	ncore_per_node=24
	nomp=12      #-c
fi

if [ $NERSC_HOST == "cori" ]; then
	nsockets=2
	ncore_per_socket=16 #two 16 haswell
	ncore_per_node=32
	nomp=16      #-c
fi