import ctypes
c_lib = ctypes.CDLL('../build/liblpmatrixgrounder.so')


def make_clist(lst):
    return (ctypes.c_char_p * len(lst))(*[x.encode() for x in lst])

filename = "../samples/example.lp".encode('utf-8')
c_lib.ground_single.argtypes = [ctypes.c_char_p]
c_lib.ground_single.restype = ctypes.c_char_p
lp = c_lib.ground_single(filename)

print(lp)

listof_files = ["../samples/model_graphcoloring.lp", "../samples/graph_6nodes.lp"]
c_lib.ground_list.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_size_t]
c_lib.ground_list.restype = ctypes.c_char_p
lp = c_lib.ground_list(make_clist(listof_files), len(listof_files))

print(lp)