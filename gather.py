

"""MPI sendrecv and gatherv example.

   Numerical DE solvers are often good candidates for parallelism
   because differential operators are local operators.  As such, the
   spatial domain of a given problem can be broken up into several
   subdomains and distributed amoungst multiple processors.  At each
   time step, the numerical solution needs to be communicated between
   processors only at the boundaries of the subdomains.

   This example shows how to break up a one-dimensional domain into
   several parts and communicate the solution at the boundaries.

   Through out this example, the root process (rank zero) is a
   marshaller and all other processes are workers.  In practise, the
   marshaller would compute the initial condition, distribute it to
   the workers, and collect the solution at the conclusion of the
   computation (or periodically during the computation) in order to
   dump it to disk.  Once the workers receive the initial condition
   they communicate amoungst themselves.

   Finally, the workers use ``ghost`` numbers of cells beyond their
   subdomain to compute the solution within their domain.  For
   example, if there is a fourth order differential operator in the
   DE, one would set ``ghost`` to 3 (or 4).

   In this example, we initialise the solution with the cell number,
   and the only computation we perform is multiplying by -1.

"""

import numpy as np
from mpi4py import MPI as mpi


######################################################################
# parameters

N     = 20                             # number of cells (grid points)
ghost = 1                              # number of ghost cells


######################################################################
# init mpi

comm    = mpi.COMM_WORLD
workers = comm.Get_size() - 1
rank    = comm.Get_rank()

if workers == 0:
    raise ValueError, 'the number of workers is zero'

#
# Initialise ``slices`` and ``gathers``:
#
# * ``slices`` - Dictionary (keyed by the rank of the worker) of index
#                arrays that describe which slice of the solution
#                (initial condition) is distributed to each worker.
#                This includes the ghost cells.
#
# * ``gathers`` - List of lengths of the above slices (excluding the
#                 ghost cells) that will be used to gather the
#                 solutions.  The first entry corresponds to the
#                 marshaller, which doesn't send anything during the
#                 gather operation.
#

slices  = {}
gathers = [ 0 ]

if workers == 1:
    slices[1] = range(N)
    gathers.append(N)

else:
    slices[1] = range(0, N/workers+ghost)
    gathers.append(N/workers)

    for i in range(2,workers):
        slices[i] = range((i-1)*(N/workers)-ghost, i*(N/workers)
+ghost)
        gathers.append(N/workers)

    slices[workers] = range(N-(N/workers + ghost + N%workers), N)
    gathers.append(N/workers + N%workers)

#
# Initialise ``right``, ``left``, ``right_worker`` and ``left_worker``
#
# * ``right`` - Index of the right-most non-ghost cell in the worker's
#               solution.
#
# * ``left`` - Index of the left-most non-ghost cell in the worker's
#              solution.
#

right = N
left  = 0

right_worker = mpi.PROC_NULL
left_worker  = mpi.PROC_NULL

if rank != 0 and workers > 1:
    M = len(slices[rank])

    if rank == 1:
        right = M - ghost
        left  = 0

        right_worker = rank + 1

    elif rank == workers:
        right = M
        left  = ghost

        left_worker = rank - 1

    else:
        right = M - ghost
        left  = ghost

        right_worker = rank + 1
        left_worker  = rank - 1


######################################################################
# allocate

if rank == 0:
    solution = np.zeros(N)                  # full domain
else:
    solution = np.zeros(len(slices[rank]))  # local subdomain


######################################################################
# compute initial condition and distribute

# initial condition (replace this)
solution = np.linspace(1.0, 1.0*N, N)

# distribute
if rank == 0:
    for i in range(1,workers+1):
        comm.Send(solution[slices[i]], dest=i)
else:
    comm.Recv(solution, source=0)


######################################################################
# compute solution

if rank != 0:

    # loop over time (replace this)
    for t in (0.0,):

        # compute solution (replace this)
        for i in range(left, right):
            solution[i] = -solution[i]

        # send non-ghost cells to right_worker's left ghost cells
        comm.Sendrecv(solution[right-ghost:right], right_worker, 0,
                      solution[:left],             left_worker,  0)

        # send non-ghost cells to left_worker's right ghost cells
        comm.Sendrecv(solution[left:left+ghost],   left_worker,  0,
                      solution[right:],            right_worker, 0)


######################################################################
# gather solution

if rank != 0:
    recv = None
    send = solution[left:right]
else:
    recv = solution
    send = None

comm.Gatherv(sendbuf=[send, mpi.DOUBLE],
             recvbuf=[recv, (gathers, None), mpi.DOUBLE],
             root=0)


######################################################################
# dump solution

if rank == 0:
    print recv














from mpi4py import MPI

#pprint = MPI._pprint
MPI.size = MPI.WORLD_SIZE
MPI.rank = MPI.WORLD_RANK

def pprint(str="", end="\n", comm=MPI.COMM_WORLD):
    """Print for MPI parallel programs: Only rank 0 prints *str*."""
    if comm.rank == 0:
        print str+end, 

# data in all processes
data = []
for i in xrange(MPI.size):
    data += [ MPI.size * 10 + MPI.rank ]


# print input data
msg = "[%d] input:  %s" % (MPI.rank, data)


pprint(msg)


# alltoall
data = MPI.WORLD.Alltoall(data)


# print result data
msg = "[%d] result: %s" % (MPI.rank, data)
pprint(msg)
