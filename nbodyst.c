#include <stdio.h>
#include <math.h>
#include <stdlib.h>


typedef struct{
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
} particle;

main () {

    int num_p = 1000;
    particle particles[num_p];
    
    FILE *fout = fopen("newsave.txt", "w");
    FILE *fin = fopen("data.dat", "rb");
    int l;
    for(l = 0; l < num_p; l++){
        fscanf(fin,"%lf %lf %lf", &particles[l].x, &particles[l].y, &particles[l].z);
        particles[l].vx = particles[l].vy = particles[l].vz =0;
        particles[l].mass = 2.0;
    }
    printf("test: %f, %f, %f\n", particles[1].x, particles[1].y, particles[1].z);
    
    float dt = 1.0;
    double G = 6.67384 * pow(10,-8);
    int i;
    int j;
    double ax, ay, az, dx, dy, dz, xi, yi, zi, r, invr, invr3, f, eps;
    eps = .01;
    int T = 100000;
    int t;

    
    for(i=0; i < num_p; i++) {
        fprintf(fout, "%f %f %f\n", particles[i].x, particles[i].y, particles[i].z);
    }
    
    for(t=0; t < T; t++){
        for(i=0; i < num_p; i++) {
            xi = particles[i].x;
            yi = particles[i].y;
            zi = particles[i].z;
            ax = 0.0;
            ay = 0.0;
            az = 0.0;
            for(j=0; j < num_p; j++) {
                if(i==j) {
                    continue;
                }
                dx = particles[j].x - xi;
                dy = particles[j].y - yi;
                dz = particles[j].z - zi;
                r = dx*dx + dy*dy + dz*dz + eps;
                invr = 1.0/sqrt(dx*dx + dy*dy + dz*dz + eps);
                invr3 = invr*invr*invr;
                f = G*particles[j].mass*invr3;
                //printf("G %f, mass %f, invr3 %f", G, particles[j].mass, invr3);
                //printf("f: %f\n", f * pow(10, 8));
                ax += f*dx;
                ay += f*dy;
                az += f*dz;
            }
            particles[i].vx += dt*ax;
            particles[i].vy += dt*ay;
            particles[i].vz += dt*az;
        }
        if(t%100 == 0) {
            printf("Timestep: %d\n", t);
        }
        for(i=0; i < num_p; i++) {
            particles[i].x += dt*particles[i].vx;
            particles[i].y += dt*particles[i].vy;
            particles[i].z += dt*particles[i].vz;
            if((t -1)%100 == 0) {
                fprintf(fout, "%f %f %f\n", particles[i].x, particles[i].y, particles[i].z);
            }
        }
        //if(t%100 == 0) {
        //    fprintf(fout, "\n");
        //}
        /*for(i=0; i < n; i++) {
            printf("Particle %d: Xpos: %f, Ypos: %f, Zpos: %f\n", i, xnew[i], ynew[i], znew[i]);
        }
        */
    }
    fclose(fout);
}
