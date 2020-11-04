#include "stepper.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>

//ldoc on
/**
 * ## Implementation
 *
 * ### Structure allocation
 */

central2d_t* central2d_init(float w, float h, int nx, int ny,
                            int block_nx, int block_ny,
			    int block_x, int block_y,
                            int nfield, flux_t flux, speed_t speed,
                            float cfl)
{
    // We extend to a four cell buffer to avoid BC comm on odd time steps
    int ng = 4;

    central2d_t* sim = (central2d_t*) malloc(sizeof(central2d_t));
    sim->nx = nx / block_nx;
    sim->ny = ny / block_nx;
    sim->block_nx = block_nx;
    sim->block_ny = block_ny;
    sim->block_x = block_x;
    sim->block_y = block_y;
    sim->ng = ng;
    sim->nfield = nfield;
    sim->dx = w / nx;
    sim->dy = h / ny;
    sim->flux = flux;
    sim->speed = speed;
    sim->cfl = cfl;

    int nx_all = sim->nx + 2*ng;
    int ny_all = sim->ny + 2*ng;
    int nc = nx_all * ny_all;
    int N  = nfield * nc;
    sim->u  = (float*) malloc((4*N + 6*nx_all)* sizeof(float));
    sim->v  = sim->u +   N;
    sim->f  = sim->u + 2*N;
    sim->g  = sim->u + 3*N;
    sim->scratch = sim->u + 4*N;

    return sim;
}


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
    return (k*ny_all+(ng+iy))*nx_all+(ng+ix);
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

static inline
void copy_subgrid2(float* restrict dst,
                   const float* restrict src,
                   int nx, int ny,
                   int dst_stride, int src_stride)
{
    for (int iy = 0; iy < ny; ++iy)
        for (int ix = 0; ix < nx; ++ix)
            dst[iy*dst_stride+ix] = src[iy*src_stride+ix];
}

void central2d_gather(central2d_t* full_sim, central2d_t* sim)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int nx = sim->nx;
    int ny = sim->ny;
    int ng = sim->ng;
    int st = nx + 2*ng;
    int nf = sim->nfield;
    int nc = (nx + 2*ng) * (ny + 2*ng);
    int N  = nf * nc;

    float* tmpu;
    if(rank == 0){
        int full_nx = full_sim->nx;
        int full_ny = full_sim->ny;
        int full_ng = full_sim->ng;
        int full_nf = full_sim->nfield;
	int full_nc = (full_nx + 2*full_ng) * (full_ny + 2*full_ng);
	int full_N  = full_nf * full_nc;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
        tmpu = (float*) malloc(N * size * sizeof(float));
    }

    MPI_Gather(sim->u, N, MPI_FLOAT, tmpu, N, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(rank == 0){
        int full_nx = full_sim->nx;
        int full_ny = full_sim->ny;
        int full_ng = full_sim->ng;
	int full_st = full_nx + 2*full_ng;
        int world;
	MPI_Comm_size(MPI_COMM_WORLD, &world);
        for(int k = 0; k < nf; k++){
            for(int i = 0; i < world; i++){
                int ibx = i % sim->block_nx;
		int iby = i / sim->block_nx;
		int ix  = ibx * sim->nx;
		int iy  = iby * sim->ny;
                float* dst = full_sim->u + central2d_offset(full_sim, k, ix, iy);
		float* src = tmpu + i*N + central2d_offset(sim, k, 0, 0);
		copy_subgrid2(dst, src, nx, ny, full_st, st);
            }
        }
    }

    if(rank == 0)
        free(tmpu);
}

void central2d_periodic(float* restrict u,
                        int nx, int ny, int ng,
                        int block_nx, int block_ny,
                        int block_x, int block_y, int nfield)
{
    int stride = nx + 2*ng;
    int N = (nx + 2*ng) * (ny + 2*ng);
    int dst_x, dst_y, dst_rank;
    int src_x, src_y, src_rank;
    float* sbuf;
    float* rbuf;

    sbuf = (float*) malloc(ng * ny * nfield * sizeof(float));
    rbuf = (float*) malloc(ng * ny * nfield * sizeof(float));

    // Send to right, receive from left
    dst_x = (block_x + 1) % block_nx;
    dst_y = block_y;
    src_x = (block_x + block_nx - 1) % block_nx;
    src_y = block_y;
    dst_rank = dst_y*block_nx + dst_x;
    src_rank = src_y*block_nx + src_x;
    for(int i = 0; i < nfield; i++){
        float* usrc = u + i*N + ng*stride + nx;
        float* sbufi = sbuf + i*(ng*ny);
        copy_subgrid2(sbufi, usrc, ng, ny, ng, stride);
    }
    MPI_Sendrecv(sbuf, nfield*(ny*ng), MPI_FLOAT, dst_rank, 0,
                 rbuf, nfield*(ny*ng), MPI_FLOAT, src_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int i = 0; i < nfield; i++){
        float* udst = u + i*N + ng*stride;
	float* rbufi = rbuf + i*(ng*ny);
        copy_subgrid2(udst, rbufi, ng, ny, stride, ng);
    }

    // Send to left, receive from right
    dst_x = (block_x + block_nx - 1) % block_nx;
    dst_y = block_y;
    src_x = (block_x + 1) % block_nx;
    src_y = block_y;
    dst_rank = dst_y*block_nx + dst_x;
    src_rank = src_y*block_nx + src_x;
    for(int i = 0; i < nfield; i++){
        float* usrc = u + i*N + ng*stride + ng;
        float* sbufi = sbuf + i*(ng*ny);
        copy_subgrid2(sbufi, usrc, ng, ny, ng, stride);
    }
    MPI_Sendrecv(sbuf, nfield*(ny*ng), MPI_FLOAT, dst_rank, 0,
                 rbuf, nfield*(ny*ng), MPI_FLOAT, src_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int i = 0; i < nfield; i++){
        float* udst = u + i*N + ng*stride + nx + ng;
	float* rbufi = rbuf + i*(ng*ny);
        copy_subgrid2(udst, rbufi, ng, ny, stride, ng);
    }

    free(sbuf);
    free(rbuf);

    sbuf = (float*) malloc(stride * ng * nfield * sizeof(float));
    rbuf = (float*) malloc(stride * ng * nfield * sizeof(float));

    // Send to bottom, receive from top
    dst_x = block_x;
    dst_y = (block_y + 1) % block_ny;
    src_x = block_x;
    src_y = (block_y + block_ny - 1) % block_ny;
    dst_rank = dst_y*block_nx + dst_x;
    src_rank = src_y*block_nx + src_x;
    for(int i = 0; i < nfield; i++){
        float* usrc = u + i*N + ny*stride;
        float* sbufi = sbuf + i*(stride*ng);
        copy_subgrid(sbufi, usrc, stride, ng, stride);
    }
    MPI_Sendrecv(sbuf, nfield*(stride*ng), MPI_FLOAT, dst_rank, 0,
                 rbuf, nfield*(stride*ng), MPI_FLOAT, src_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int i = 0; i < nfield; i++){
        float* udst = u + i*N;
        float* rbufi = rbuf + i*(stride*ng);
        copy_subgrid(udst, rbufi, stride, ng, stride);
    }

    // Send to top, receive from bottom
    dst_x = block_x;
    dst_y = (block_y + block_ny - 1) % block_ny;
    src_x = block_x;
    src_y = (block_y + 1) % block_ny;
    dst_rank = dst_y*block_nx + dst_x;
    src_rank = src_y*block_nx + src_x;
    for(int i = 0; i < nfield; i++){
        float* usrc = u + i*N + ng*stride;
        float* sbufi = sbuf + i*(stride*ng);
        copy_subgrid(sbufi, usrc, stride, ng, stride);
    }
    MPI_Sendrecv(sbuf, nfield*(stride*ng), MPI_FLOAT, dst_rank, 0,
                 rbuf, nfield*(stride*ng), MPI_FLOAT, src_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for(int i = 0; i < nfield; i++){
        float* udst = u + i*N + (ny + ng)*stride;
        float* rbufi = rbuf + i*(stride*ng);
        copy_subgrid(udst, rbufi, stride, ng, stride);
    }

    free(sbuf);
    free(rbuf);
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
                   int block_nx, int block_ny,
                   int block_x, int block_y,
                   int nfield, flux_t flux, speed_t speed,
                   float tfinal, float dx, float dy, float cfl)
{
    int nstep = 0;
    int nx_all = nx + 2*ng;
    int ny_all = ny + 2*ng;
    bool done = false;
    float t = 0;
    while (!done) {
        float cxy[2] = {1.0e-15f, 1.0e-15f};
        central2d_periodic(u, nx, ny, ng, block_nx, block_ny,
                           block_x, block_y, nfield);
        speed(cxy, u, nx_all * ny_all, nx_all * ny_all);
	//printf("cxy: %f %f\n", cxy[0], cxy[1]);
        float my_dt = cfl / fmaxf(cxy[0]/dx, cxy[1]/dy);
	float dt;
	MPI_Allreduce(&my_dt, &dt, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
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
    return central2d_xrun(sim->u, sim->v, sim->scratch,
                          sim->f, sim->g,
                          sim->nx, sim->ny, sim->ng,
                          sim->block_nx, sim->block_ny,
                          sim->block_x, sim->block_y,
                          sim->nfield, sim->flux, sim->speed,
                          tfinal, sim->dx, sim->dy, sim->cfl);
}
