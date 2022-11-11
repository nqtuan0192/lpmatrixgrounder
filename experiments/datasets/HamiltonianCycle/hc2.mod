using CP;

{string} node = ... ;
tuple edge_type {string edge_0;string edge_1;};
{edge_type}edge = ... ;
{string} start = ... ;
int reached_start[node] = ... ;
int hc[edge] = ... ;
int reached[node] = ... ;
subject to{
	forall(x0 in node){
		(1)>0 => (reached_start[x0])>0;
	}
};

