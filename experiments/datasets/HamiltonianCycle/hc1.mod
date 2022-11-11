using CP;

{string} node = ... ;
tuple edge_type {string edge_0;string edge_1;};
{edge_type}edge = ... ;
{string} start = ... ;
int reached_start[node];
dvar boolean hc[edge];
dvar boolean reached[node];

subject to{
	forall(<x0,x2> in edge, <x0,x1> in edge: x1!=x2){
		((hc[<x0,x1>])*(hc[<x0,x2>])*1)>0 => false;
	}
	forall(<x0,x1> in edge, <x2,x1> in edge: x0!=x2){
		((hc[<x0,x1>])*(hc[<x2,x1>])*1)>0 => false;
	}
	forall(x0 in node){
		(1)>0 => (reached[x0])>0;
	}
	forall(x0 in node){
		(reached[x0])>0 => (sum(<y0,x0> in edge) ((hc[<y0,x0>])))>0;
		(sum(<y0,x0> in edge) ((hc[<y0,x0>])))>0 => (reached[x0])>0;
	}
};

execute{
	for(var z0 in node)
		reached_start[z0]=0;
	for (var x0 in start){
		for (var x1 in edge){
			if(x0==x1.edge_0 && hc[x1]==1){
				reached_start[x1.edge_1]=1;
			}
		}
	}
	var modified0=true;
	while(modified0){
		modified0=false;
		for (var x2 in edge){
			if(hc[x2]==1 && reached_start[x2.edge_0]==1 && reached_start[x2.edge_1]==0){
				reached_start[x2.edge_1]=1;
				modified0=true;
			}
		}
	}
}

