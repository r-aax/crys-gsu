from split import split
from mpi4py import MPI
import numpy as np
import multiprocessing as mp
import sys
from remeshing import TongRemesher
import functools
from gsu.zones_adjacency_matrix import ZonesAdjacencyMatrix
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
import swim
np.set_printoptions(threshold=sys.maxsize)
sys.setrecursionlimit(100000)

def get_data(face, name):
    """Get data from Face.

    Parameters
    ----------
      face : Face
        face
      name : string
        name of the attribute
    
    Returns
    -------
      float
    """
    return face.Data[name]

def get_mpi_borders(rank, zam, z, g):
    """Compose dictionary of type:
        { 
          rank A: [edge1 id, edge2 id, ..., edgeN id],
          rank B: [...] 
        }
        where A and B are the ranks of neigbhoring zones.
    
    Parameters
    ----------
      rank : int
        rank of the current process
      zam : ZonesAdjacencyMatrix obj
        matrix of zones adjacency
      z : Zone
        current zone
      g : Grid
        grid
    
    Returns
    -------
      dict
    """
    line = zam.zone_cross_edges_array(rank)
    mpi_borders = {}
    loc_edges_ids = [-1] * len(g.Edges)
    for i, edge in enumerate(z.Edges):
        loc_edges_ids[edge.GloId] = i

    for (li, le) in enumerate(line):
        if le > 0:
            lz = g.Zones[li]
            mpi_borders[li] = list()
            for e in z.Edges:
                if e.is_connect_zones(z, lz):
                    mpi_borders[li].append(loc_edges_ids[e.GloId])
    return mpi_borders

def exchange(mpi_borders, send_sb, recv_sb, data):
    """Exchange data between all zones.

    Parameters
    ----------
      mpi_borders : dict
        information about connectivity between adjacent zones.
      send_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
      recv_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
      data : string
        name of attribute
    """
    put_to_superbuffer(mpi_borders, send_sb, data)
    start_exchanges(mpi_borders, send_sb, recv_sb)
    get_from_superbuffer(mpi_borders, recv_sb, data)

def put_to_superbuffer(mpi_borders, send_sb, data):
    """Put data into send superbuffer.

    Parameters
    ----------
      mpi_borders : dict
        information about connectivity between adjacent zones.
      send_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
      data : string
        name of attribute
    """
    for i, v in enumerate(mpi_borders.values()):
        for j, eid in enumerate(v):
            send_sb[i, j] = z.get_real_face(eid)[data]

def start_exchanges(mpi_borders, send_sb, recv_sb):
    """Use non-blocking MPI calls to exchage data between superbuffers.

    Parameters
    ----------
      mpi_borders : dict
        information about connectivity between adjacent zones.
      send_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
      recv_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
    """
    requests = list()
    for i in range(len(mpi_borders)):
        dest = list(mpi_borders.keys())[i]
        count = len(mpi_borders[dest])
        comm.Isend(send_sb[i][:count], dest=dest)

    for i in range(len(mpi_borders)):
        source = list(mpi_borders.keys())[i]
        count = len(mpi_borders[source])
        req = comm.Irecv(recv_sb[i][:count], source=source)
        requests.append(req)

    MPI.Request.Waitall(requests)
    comm.Ibarrier()

def get_from_superbuffer(mpi_borders, recv_sb, data):
    """Get data from recv superbuffer.

    Parameters
    ----------
      mpi_borders : dict
        information about connectivity between adjacent zones.
      recv_sb : np n dim array
        np array of shape [N_NEIGBHORS, MAX_EDGES]
      data : string
        name of attribute
    """
    for i, v in enumerate(mpi_borders.values()):
        for j, eid in enumerate(v):
            z.get_ghost_face(eid)[data] = recv_sb[i, j]

def remesh(rank, his):
    """Apply remeshing based on ice heights collected from all zones with MPI Gather.
    Store grid into file.

    Parameters
    ----------
      rank : int
        rank of the process
      his : np array
        np array of shape [N_FACES]
    """
    sendbuf = his
    sendcounts = np.array(comm.gather(sendbuf.size, root))
    print(sendcounts)
    if rank == root:
        print("sendcounts: {}, total: {}".format(sendcounts, sum(sendcounts)))
        recvbuf = np.empty(sum(sendcounts), dtype=np.float64)
    else:
        recvbuf = None

    comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)
    if rank == root:
        print("Gathered array: {}".format(recvbuf.size))

    if rank == 0:
        for i, zone in enumerate(g.Zones):
            for f, new_h in zip(zone.Faces, recvbuf[i * sendcounts[i] : (i + 1) * sendcounts[i]]):
                f.Data['Hi'] = new_h

        TongRemesher(mesh_filename=FILENAME, json_file=JSON, outdir=OUTDIR)(g)
        g.store(OUT_REMESHED_MESH)


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    root = 0

    FILENAME = '../grids/bunny.dat'
    JSON = '../grids/bunny.json'
    OUTDIR = '../'
    OUT_DAT_FILE = 'res.dat'
    CRYDIR = '.'
    OUT_REMESHED_MESH = OUTDIR + 'remeshed.dat'

    g = split(FILENAME, cry_dir=CRYDIR, out_dat_file=OUTDIR + OUT_DAT_FILE, split_strategy='n{}'.format(size))
    zam = ZonesAdjacencyMatrix(g.Edges, g.Zones)
    z = g.Zones[rank]

    pool = mp.Pool(processes=6)
    HIS = np.array([pool.map(functools.partial(get_data, name='Hi'), z.Faces)], dtype=np.float64)
    T = np.array([pool.map(functools.partial(get_data, name='T'), z.Faces)], dtype=np.float64)
    HW = np.array([pool.map(functools.partial(get_data, name='Hw'), z.Faces)], dtype=np.float64)

    swim.hi(T[0], HIS[0], HW[0])

    if size > 1:
        mpi_borders = get_mpi_borders(rank, zam, z, g)
        max_length = max(list(map(len, mpi_borders.values())))

        # todo numpy не позволяет делать массив с массивами разной длины, поэтому буфферы строятся по максимальной границе.
        send_superbuffer = np.empty((len(mpi_borders), max_length))
        recv_superbuffer = np.empty((len(mpi_borders), max_length))

        for f in z.Faces:
            f.Zone = z.Id

        exchange(mpi_borders, send_superbuffer, recv_superbuffer, 'T')
        exchange(mpi_borders, send_superbuffer, recv_superbuffer, 'Hi')
        exchange(mpi_borders, send_superbuffer, recv_superbuffer, 'Hw')

    remesh(rank, HIS)
