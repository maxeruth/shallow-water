#include "stepper.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

//ldoc on
/**
 * ## Implementation
 *
 * ### Structure allocation
 */

central2d_t* central2d_init(float w, float h, int nx_total, int ny_total,
                            int ng, int NX, int NY,
                            int nfield, flux_t flux, speed_t speed,
                            float cfl)
{
	// For definitions of the values, see stepper.h
    central2d_t* sim = (central2d_t*) malloc(sizeof(central2d_t));
    
    // Get number of processors, make sure that NX*NY is the world size
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (NX*NY != world_size){
        fprintf(stderr, "Number of processors in grid (%d) does not match world size (%d)\n", NX*NY, world_size);
        exit(-1);
    }
    
    sim->world_size = world_size;
		
    // Get rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    sim->rank = rank;
    
    // Get the positions (rank = Y * NY + X)
    int X,Y;
    Y = rank/NY; X = rank%NY;
    
    // Get the offsets for loops
    int x0,y0;
    x0 = X * (nx_total/NX);
    y0 = Y * (ny_total/NY);
    
    sim->x0 = x0;
    sim->y0 = y0;
    
    // Get the total number of points "owned" in each direction
    int nx,ny;
    nx = ((nx_total/NX) < (nx_total-x0)) ? (nx_total/NX) : (nx_total-x0);
    ny = ((ny_total/NY) < (ny_total-y0)) ? (ny_total/NY) : (ny_total-y0);
    
    sim->nx = nx;
    sim->ny = ny;
    
    sim->NX = NX;
    sim->NY = NY;
    sim->ng = ng;
    sim->nfield = nfield;
    sim->dx = w/nx_total;
    sim->dy = h/ny_total;
    sim->flux = flux;
    sim->speed = speed;
    sim->cfl = cfl;
    
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    int nc = nx_all * ny_all;
    int N  = nfield * nc;
    sim->u  = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
    sim->v  = sim->u +   N;
    sim->f  = sim->u + 2*N;
    sim->g  = sim->u + 3*N;
    sim->scratch = sim->u + 4*N;
    
    sim->neighbors[0] = Y*NX + (X+NX-1)%NX;        // left
    sim->neighbors[1] = Y*NX + (X+NX+1)%NX;        // right
    sim->neighbors[2] = ((Y+1)%NY)*NX + X;         // top
    sim->neighbors[3] = ((Y+NY-1)%NY)*NX + X;      // bottom

	printf("\nsim set up\n"
	       //"rank = %d, nfield = %d, nx = %d, ny = %d, ng = %d,\n"
	       //"rank = %d, dx = %g, dy = %g, cfl = %g, world_size = %d,\n"
	       //"rank = %d, NX = %d, NY = %d, x0 = %d, y0 = %d,\n"
	       "rank = %d, bottom_neighbor = %d, top_neighbor = %d,\n"
	       "rank = %d, left_neighbor = %d, right_neighbor = %d\n\n",
	       //sim->rank,sim->nfield,sim->nx,sim->ny,sim->ng,
	       //sim->rank,sim->dx,sim->dy,sim->cfl,sim->world_size,
	       //sim->rank,sim->NX,sim->NY,sim->x0,sim->y0,
	       sim->rank,sim->neighbors[3],sim->neighbors[2],
	       sim->rank,sim->neighbors[0],sim->neighbors[1]);
	
	// printf("Just to make sure: rank - Y*NX+X = %d\n",rank - (Y*NX+X));
	
    return sim;
}


void copy_basic_info(int nx, int ny, central2d_t* sim, central2d_t* full_sim){
	full_sim->nfield = sim->nfield;
	full_sim->nx = nx;
	full_sim->ny = ny;
	full_sim->ng = sim->ng;
	full_sim->dx = sim->dx;
	full_sim->dy = sim->dy;
	full_sim->cfl = sim->cfl;
	full_sim->rank = sim->rank;
	full_sim->world_size = sim->world_size; 
	full_sim->NX = sim->NX;
	full_sim->NY = sim->NY;
	full_sim->x0 = 0;
	full_sim->y0 = 0;
    for (int i = 0; i < 4; ++i){
        full_sim->neighbors[i] = 0;
    }
	
	full_sim->flux = sim->flux;
	full_sim->speed = sim->speed;
	
	int nx_all = nx + 2*full_sim->ng;
    int ny_all = ny + 2*full_sim->ng;
    int nc = nx_all * ny_all;
    int N  = full_sim->nfield * nc;
	full_sim->u  = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
	full_sim->v  = sim->u +   N;
    full_sim->f  = sim->u + 2*N;
    full_sim->g  = sim->u + 3*N;
    full_sim->scratch = sim->u + 4*N;
	
}

// Copy the data from the source into the full_sim
void copy_u(float* u, int source_nx, int source_ny, 
            int source_x0, int source_y0, central2d_t* full_sim){
    printf("Entering copy_u w/ nx = %d, ny = %d, x0 = %d, y0 = %d\n",source_nx,source_ny,source_x0,source_y0);
    int ng = full_sim->ng;
    int nx_all = source_nx + 2*ng;
    int ny_all = source_ny + 2*ng;
    int N = nx_all * ny_all;

    for (int iy = 0; iy < source_ny; ++iy){
        for (int ix = 0; ix < source_nx; ++ix){
            int iu = (ng+iy)*nx_all + (ng+ix);
            full_sim->u[central2d_offset(full_sim,0,source_x0 + ix,source_y0 + iy)] = u[iu];
            full_sim->u[central2d_offset(full_sim,1,source_x0 + ix,source_y0 + iy)] = u[N + iu];
            full_sim->u[central2d_offset(full_sim,2,source_x0 + ix,source_y0 + iy)] = u[2*N + iu];
            full_sim->u[central2d_offset(full_sim,3,source_x0 + ix,source_y0 + iy)] = u[3*N + iu];
            /*printf("x0-%d y0-%d snx-%d, sny-%d copy_u: ix = %d, iy = %d, iu = %d, iu_full = %d u = %g\n",
                   source_x0,source_y0,source_nx,source_ny,ix,iy,iu,central2d_offset(full_sim,0,source_x0+ix,source_y0+iy), 
                   full_sim->u[central2d_offset(full_sim,0,source_x0+ix,source_y0+iy)]);*/
        }
    }
}


// Send data from sim->u to a destination (probably the rank 0 node)
void send_full_u(int destination, central2d_t* sim){
    int nx_all = sim->nx + 2*sim->ng;
    int ny_all = sim->ny + 2*sim->ng;
    int nc = nx_all * ny_all;
    int N  = sim->nfield * nc;
    printf("Sending u from %d to %d\n",sim->rank,destination);
    MPI_Send(sim->u, 4*N + 6*nx_all, MPI_FLOAT, destination, 0, MPI_COMM_WORLD);
}

// Receive data from send_full_u. Assumes full_sim is the full simulation for writes and solution checking
void recv_full_u(int source, central2d_t* full_sim){
    // Calculate the amount of data that must be coming from source:
    int X,Y;
    int NX,NY;
    NX = full_sim->NX; 
    NY = full_sim->NY;
    Y = source/NY; X = source%NY;
    
    int x0,y0;
    x0 = X * (full_sim->nx/NX);
    y0 = Y * (full_sim->ny/NY);
    
    int nx,ny;
    nx = ( ( full_sim->nx / NX ) < ( full_sim->nx - x0 ) ) ? ( full_sim->nx / NX ) : ( full_sim->nx - x0 );
    ny = ( ( full_sim->ny / NY ) < ( full_sim->ny - y0 ) ) ? ( full_sim->ny / NY ) : ( full_sim->ny - y0 );
    
    int ng = full_sim->ng;
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    int nc = nx_all * ny_all;
    int N  = full_sim->nfield * nc;
    
    // Create a temporary place to receive all of the data
    float* tmpu = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
    
    // Receive the data from the source node
    printf("Receiving from %d to %d\n",source,full_sim->rank);
    MPI_Recv(tmpu, 4*N + 6*nx_all, MPI_FLOAT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
    // Copy the data into the full solution array
    copy_u(tmpu, nx, ny, x0, y0, full_sim);
    
    free(tmpu);
}


void gather_sol(central2d_t* sim, central2d_t* full_sim){
	printf("r%d Entering gather_sol\n",sim->rank);
	if(sim->rank == 0){
		// Copy sim into full_sim
		copy_u(sim->u, sim->nx, sim->ny, sim->x0, sim->y0, full_sim);
		
		for(int ii = 1; ii < sim->world_size; ++ii){
			recv_full_u(ii,full_sim); // Receive from all other nodes, synthesize into 
		}
		
	}
	else{
		send_full_u(0,sim); // Send to node 0	
	}
}

/*
central2d_t* copy_subdomain(central2d_t* sim, int num_domain)
{
    //int BLOCK_X = 4;
    //int BLOCK_Y = 4;
    central2d_t* sim_sub = (central2d_t*) malloc(sizeof(central2d_t));
    sim_sub->nx = sim->nx/BLOCK_X; // BLOCK_X number of blocks in x direction
    sim_sub->ny = sim->ny/BLOCK_Y; // BLOCK_Y number of blocks in y direction
    sim_sub->ng = sim->ng;
    sim_sub->nfield = sim->nfield;
    sim_sub->dx = sim->dx;
    sim_sub->dy = sim->dy;
    sim_sub->flux = sim->flux;
    sim_sub->speed = sim->speed;
    sim_sub->cfl = sim->cfl;

    int N = sim->nx + 2*sim->ng; // size of original matrix plus ghost cells
    int M = sim_sub->nx + 2*sim_sub->ng; // size of blocking matrix plus ghost cells
    int nc_block = M * M;
    int NN  = sim_sub->nfield * nc_block;
    sim_sub->u  = (float*) malloc((4*NN + 6*M)* sizeof(float));

    // coordinates of the blocking matrix in original matrix without ghost cells
    int bx = num_domain/BLOCK_X;
    int by = num_domain%BLOCK_X;
    int i_bx = bx * sim_sub->nx;
    int i_by = by * sim_sub->ny;

    // copy params to subdomain for parallelization
    for (int iy = 0; iy < M; ++iy){
        for (int ix = 0; ix < M; ++ix){
            sim_sub->u[iy*M+ix] = (sim->u + central2d_offset(sim, 0, i_bx, i_by) - N)[iy*N + ix];
            sim_sub->u[iy*M+ix + NN] = (sim->u + central2d_offset(sim, 1, i_bx, i_by) - N)[iy*N + ix];
            sim_sub->u[iy*M+ix + 2*NN] = (sim->u + central2d_offset(sim, 2, i_bx, i_by) - N)[iy*N + ix];
            sim_sub->u[iy*M+ix + 3*NN] = (sim->u + central2d_offset(sim, 3, i_bx, i_by) - N)[iy*N + ix];
        }
    }

    sim_sub->v  = sim_sub->u +   NN;
    sim_sub->f  = sim_sub->u + 2*NN;
    sim_sub->g  = sim_sub->u + 3*NN;
    sim_sub->scratch = sim_sub->u + 4*NN;

    return sim_sub;
} */



void central2d_free(central2d_t* sim)
{
    free(sim->u);
    free(sim);
}


int central2d_offset(central2d_t* sim, int k, int ix, int iy)
{
    int nx = sim->nx, ny = sim->ny, ng = sim->ng;
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    return (k*ny_all + (ng+iy))*nx_all + (ng+ix);
}


/**
 * ### Boundary conditions
 *
 * In finite volume methods, boundary conditions are typically applied by
 * setting appropriate values in ghost cells.  For our framework, we will
 * apply periodic boundary conditions; that is, waves that exit one side
 * of the domain will enter from the other side.
 *
 * We apply the conditions by assuming that the cells with coordinates
 * `nghost <= ix <= nx+nghost` and `nghost <= iy <= ny+nghost` are
 * "canonical", and setting the values for all other cells `(ix,iy)`
 * to the corresponding canonical values `(ix+p*nx,iy+q*ny)` for some
 * integers `p` and `q`.
 */

static inline
void copy_subgrid(float* restrict dst,
                  const float* restrict src,
                  int nx, int ny, int stride)
{
    for (int iy = 0; iy < ny; ++iy)
        for (int ix = 0; ix < nx; ++ix)
            dst[iy*stride+ix] = src[iy*stride+ix];
}

void central2d_periodic(float* restrict u,
                        int nx, int ny, int ng, int nfield)
{
    // Stride and number per field
    int s = nx + 2*ng;
    int field_stride = (ny+2*ng)*s;

    // Offsets of left, right, top, and bottom data blocks and ghost blocks
    int l = nx,   lg = 0;
    int r = ng,   rg = nx+ng;
    int b = ny*s, bg = 0;
    int t = ng*s, tg = (ny+ng)*s;

    // Copy data into ghost cells on each side
    for (int k = 0; k < nfield; ++k) {
        float* uk = u + k*field_stride;
        copy_subgrid(uk+lg, uk+l, ng, ny+2*ng, s);
        copy_subgrid(uk+rg, uk+r, ng, ny+2*ng, s);
        copy_subgrid(uk+tg, uk+t, nx+2*ng, ng, s);
        copy_subgrid(uk+bg, uk+b, nx+2*ng, ng, s);
    }
}

void central2d_periodic_send(float* restrict u,
                             int nx, int ny, int ng, int nfield,
                             int rank, int neighbors[4])
{
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    int nc = nx_all * ny_all;
    int N  = nfield * nc;
    
    for (int i = 0; i < 4; ++i){
        printf("Sending u from %d to %d\n",rank, neighbors[i]);
        MPI_Send(u, 4*N + 6*nx_all, MPI_FLOAT, neighbors[i], i, MPI_COMM_WORLD);
    }
}

void central2d_periodic_recv(float* restrict u,
                             int nx, int ny, int ng, int nfield,
                             int rank, int neighbors[4])
{
    // Stride and number per field
    int s = nx + 2*ng;
    int field_stride = (ny+2*ng)*s;
    int N = nfield * field_stride;

    // Offsets of left, right, top, and bottom data blocks and ghost blocks
    int block_offset[4] = {nx, ng, ng*s, ny*s};
    int gblock_offset[4] = {0, nx+ng, (ny+ng)*s, 0};
    int nxk[4] = {ng, ng, nx+2*ng, nx+2*ng};
    int nyk[4] = {ny+2*ng, ny+2*ng, ng, ng};
    
    for (int i = 0; i < 4; i++){
    
    // Create a temporary place to receive all of the data
    float* tmpu = (float*) malloc((4*N + 6*s)* sizeof(float));
    
    // Receive the data from the source node
    printf("Receiving from %d to %d\n",neighbors[i], rank);
    MPI_Recv(tmpu, 4*N + 6*s, MPI_FLOAT, neighbors[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // copy data from recv'd data to ghost blocks
    for (int k = 0; k < nfield; ++k) {
        float* uk = u + k*field_stride;
        float* recv_uk = tmpu + k*field_stride;
        copy_subgrid(uk+gblock_offset[i], recv_uk+block_offset[i], nxk[i], nyk[i], s);
        free(recv_uk);
    }
    free(tmpu);
    }
}


/**
 * ### Derivatives with limiters
 *
 * In order to advance the time step, we also need to estimate
 * derivatives of the fluxes and the solution values at each cell.
 * In order to maintain stability, we apply a limiter here.
 *
 * The minmod limiter *looks* like it should be expensive to computer,
 * since superficially it seems to require a number of branches.
 * We do something a little tricky, getting rid of the condition
 * on the sign of the arguments using the `copysign` instruction.
 * If the compiler does the "right" thing with `max` and `min`
 * for floating point arguments (translating them to branch-free
 * intrinsic operations), this implementation should be relatively fast.
 */


// Branch-free computation of minmod of two numbers times 2s
static inline
float xmin2s(float s, float a, float b) {
    float sa = copysignf(s, a);
    float sb = copysignf(s, b);
    float abs_a = fabsf(a);
    float abs_b = fabsf(b);
    float min_abs = (abs_a < abs_b ? abs_a : abs_b);
    return (sa+sb) * min_abs;
}


// Limited combined slope estimate
static inline
float limdiff(float um, float u0, float up) {
    const float theta = 2.0;
    const float quarter = 0.25;
    float du1 = u0-um;   // Difference to left
    float du2 = up-u0;   // Difference to right
    float duc = up-um;   // Twice centered difference
    return xmin2s( quarter, xmin2s(theta, du1, du2), duc );
}


// Compute limited derivs
static inline
void limited_deriv1(float* restrict du,
                    const float* restrict u,
                    int ncell)
{
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i-1], u[i], u[i+1]);
}


// Compute limited derivs across stride
static inline
void limited_derivk(float* restrict du,
                    const float* restrict u,
                    int ncell, int stride)
{
    assert(stride > 0);
    for (int i = 0; i < ncell; ++i)
        du[i] = limdiff(u[i-stride], u[i], u[i+stride]);
}


/**
 * ### Advancing a time step
 *
 * Take one step of the numerical scheme.  This consists of two pieces:
 * a first-order corrector computed at a half time step, which is used
 * to obtain new $F$ and $G$ values; and a corrector step that computes
 * the solution at the full step.  For full details, we refer to the
 * [Jiang and Tadmor paper][jt].
 *
 * The `compute_step` function takes two arguments: the `io` flag
 * which is the time step modulo 2 (0 if even, 1 if odd); and the `dt`
 * flag, which actually determines the time step length.  We need
 * to know the even-vs-odd distinction because the Jiang-Tadmor
 * scheme alternates between a primary grid (on even steps) and a
 * staggered grid (on odd steps).  This means that the data at $(i,j)$
 * in an even step and the data at $(i,j)$ in an odd step represent
 * values at different locations in space, offset by half a space step
 * in each direction.  Every other step, we shift things back by one
 * mesh cell in each direction, essentially resetting to the primary
 * indexing scheme.
 *
 * We're slightly tricky in the corrector in that we write
 * $$
 *   v(i,j) = (s(i+1,j) + s(i,j)) - (d(i+1,j)-d(i,j))
 * $$
 * where $s(i,j)$ comprises the $u$ and $x$-derivative terms in the
 * update formula, and $d(i,j)$ the $y$-derivative terms.  This cuts
 * the arithmetic cost a little (not that it's that big to start).
 * It also makes it more obvious that we only need four rows worth
 * of scratch space.
 */


// Predictor half-step
static
void central2d_predict(float* restrict v,
                       float* restrict scratch,
                       const float* restrict u,
                       const float* restrict f,
                       const float* restrict g,
                       float dtcdx2, float dtcdy2,
                       int nx, int ny, int nfield)
{
    float* restrict fx = scratch;
    float* restrict gy = scratch+nx;
    for (int k = 0; k < nfield; ++k) {
        for (int iy = 1; iy < ny-1; ++iy) {
            int offset = (k*ny+iy)*nx+1;
            limited_deriv1(fx+1, f+offset, nx-2);
            limited_derivk(gy+1, g+offset, nx-2, nx);
            for (int ix = 1; ix < nx-1; ++ix) {
                int offset = (k*ny+iy)*nx+ix;
                v[offset] = u[offset] - dtcdx2 * fx[ix] - dtcdy2 * gy[ix];
            }
        }
    }
}


// Corrector
static
void central2d_correct_sd(float* restrict s,
                          float* restrict d,
                          const float* restrict ux,
                          const float* restrict uy,
                          const float* restrict u,
                          const float* restrict f,
                          const float* restrict g,
                          float dtcdx2, float dtcdy2,
                          int xlo, int xhi)
{
    for (int ix = xlo; ix < xhi; ++ix)
        s[ix] =
            0.2500f * (u [ix] + u [ix+1]) +
            0.0625f * (ux[ix] - ux[ix+1]) +
            dtcdx2  * (f [ix] - f [ix+1]);
    for (int ix = xlo; ix < xhi; ++ix)
        d[ix] =
            0.0625f * (uy[ix] + uy[ix+1]) +
            dtcdy2  * (g [ix] + g [ix+1]);
}


// Corrector
static
void central2d_correct(float* restrict v,
                       float* restrict scratch,
                       const float* restrict u,
                       const float* restrict f,
                       const float* restrict g,
                       float dtcdx2, float dtcdy2,
                       int xlo, int xhi, int ylo, int yhi,
                       int nx, int ny, int nfield)
{
    assert(0 <= xlo && xlo < xhi && xhi <= nx);
    assert(0 <= ylo && ylo < yhi && yhi <= ny);

    float* restrict ux = scratch;
    float* restrict uy = scratch +   nx;
    float* restrict s0 = scratch + 2*nx;
    float* restrict d0 = scratch + 3*nx;
    float* restrict s1 = scratch + 4*nx;
    float* restrict d1 = scratch + 5*nx;

    for (int k = 0; k < nfield; ++k) {

        float*       restrict vk = v + k*ny*nx;
        const float* restrict uk = u + k*ny*nx;
        const float* restrict fk = f + k*ny*nx;
        const float* restrict gk = g + k*ny*nx;

        limited_deriv1(ux+1, uk+ylo*nx+1, nx-2);
        limited_derivk(uy+1, uk+ylo*nx+1, nx-2, nx);
        central2d_correct_sd(s1, d1, ux, uy,
                             uk + ylo*nx, fk + ylo*nx, gk + ylo*nx,
                             dtcdx2, dtcdy2, xlo, xhi);

        for (int iy = ylo; iy < yhi; ++iy) {

            float* tmp;
            tmp = s0; s0 = s1; s1 = tmp;
            tmp = d0; d0 = d1; d1 = tmp;

            limited_deriv1(ux+1, uk+(iy+1)*nx+1, nx-2);
            limited_derivk(uy+1, uk+(iy+1)*nx+1, nx-2, nx);
            central2d_correct_sd(s1, d1, ux, uy,
                                 uk + (iy+1)*nx, fk + (iy+1)*nx, gk + (iy+1)*nx,
                                 dtcdx2, dtcdy2, xlo, xhi);

            for (int ix = xlo; ix < xhi; ++ix)
                vk[iy*nx+ix] = (s1[ix]+s0[ix])-(d1[ix]-d0[ix]);
        }
    }
}


static
void central2d_step(float* restrict u, float* restrict v,
                    float* restrict scratch,
                    float* restrict f,
                    float* restrict g,
                    int io, int nx, int ny, int ng,
                    int nfield, flux_t flux, speed_t speed,
                    float dt, float dx, float dy)
{
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;

    float dtcdx2 = 0.5 * dt / dx;
    float dtcdy2 = 0.5 * dt / dy;

    flux(f, g, u, nx_all * ny_all, nx_all * ny_all);

    central2d_predict(v, scratch, u, f, g, dtcdx2, dtcdy2,
                      nx_all, ny_all, nfield);

    // Flux values of f and g at half step
    for (int iy = 1; iy < ny_all-1; ++iy) {
        int jj = iy*nx_all+1;
        flux(f+jj, g+jj, v+jj, nx_all-2, nx_all * ny_all);
    }

    central2d_correct(v+io*(nx_all+1), scratch, u, f, g, dtcdx2, dtcdy2,
                      ng-io, nx+ng-io,
                      ng-io, ny+ng-io,
                      nx_all, ny_all, nfield);
}


/**
 * ### Advance a fixed time
 *
 * The `run` method advances from time 0 (initial conditions) to time
 * `tfinal`.  Note that `run` can be called repeatedly; for example,
 * we might want to advance for a period of time, write out a picture,
 * advance more, and write another picture.  In this sense, `tfinal`
 * should be interpreted as an offset from the time represented by
 * the simulator at the start of the call, rather than as an absolute time.
 *
 * We always take an even number of steps so that the solution
 * at the end lives on the main grid instead of the staggered grid.
 */

static
int central2d_xrun(float* restrict u, float* restrict v,
                   float* restrict scratch,
                   float* restrict f,
                   float* restrict g,
                   int nx, int ny, int ng,
                   int nfield, flux_t flux, speed_t speed,
                   float tfinal, float dx, float dy, float cfl,
                   int rank, int neighbors[4])
{
    int nstep = 0;
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    bool done = false;
    float t = 0;
    while (!done) {
        float cxy[2] = {1.0e-15f, 1.0e-15f};
        // change central2d_periodic into data sharing from left, right, top, bottom processors
        // central2d_periodic(u, nx, ny, ng, nfield);
        printf("Data sharing start...");
        central2d_periodic_send(u, nx, ny, ng, nfield, rank, neighbors);
        central2d_periodic_recv(u, nx, ny, ng, nfield, rank, neighbors);
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Data sharing received...");
        speed(cxy, u, nx_all * ny_all, nx_all * ny_all);
        float dt = cfl / fmaxf(cxy[0]/dx, cxy[1]/dy);
        if (t + 2*dt >= tfinal) {
            dt = (tfinal-t)/2;
            done = true;
        }
        central2d_step(u, v, scratch, f, g,
                       0, nx+4, ny+4, ng-2,
                       nfield, flux, speed,
                       dt, dx, dy);
        central2d_step(v, u, scratch, f, g,
                       1, nx, ny, ng,
                       nfield, flux, speed,
                       dt, dx, dy);
        t += 2*dt;
        nstep += 2;
    }
    return nstep;
}




int central2d_run(central2d_t* sim, float tfinal)
{
    /*
    int num_domain = 16;
    for (int i = 0; i < num_domain; ++i){
    central2d_t* sim_sub = copy_subdomain(sim, i); // blocking, copy the original domain into a few sub-domains
    
    central2d_xrun(sim_sub->u, sim_sub->v, sim_sub->scratch,
                   sim_sub->f, sim_sub->g,
                   sim_sub->nx, sim_sub->ny, sim_sub->ng,
                   sim_sub->nfield, sim_sub->flux, sim_sub->speed,
                   tfinal, sim_sub->dx, sim_sub->dy, sim_sub->cfl);
    }*/
    
    return central2d_xrun(sim->u, sim->v, sim->scratch,
                          sim->f, sim->g,
                          sim->nx, sim->ny, sim->ng,
                          sim->nfield, sim->flux, sim->speed,
                          tfinal, sim->dx, sim->dy, sim->cfl,
                          sim->rank, sim->neighbors);
}
