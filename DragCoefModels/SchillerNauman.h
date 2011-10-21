Thread* thread_l = THREAD_SUB_THREAD(mix_thread, s_col); /*liquid phase - primary*/
Thread* thread_s = THREAD_SUB_THREAD(mix_thread, f_col); /*solid phase - secondary*/

real slip_x = C_U(cell, thread_l) - C_U(cell, thread_s);
real slip_y = C_V(cell, thread_l) - C_V(cell, thread_s);
real slip_z = C_W(cell, thread_l) - C_W(cell, thread_s);

/*Euclidean norm of slip velocity*/
real slip = sqrt(slip_x*slip_x + slip_y*slip_y + slip_z*slip_z);

/*relative Reynolds number for the primary phase*/
real reyp = RHO_L*slip*DIAMETER/MU_L;

/* compute particle relaxation time */
real taup = RHO_S*DIAMETER*DIAMETER/18./MU_L;

/*drag coefficient*/
real cd0;

/*Schiller-Nauman*/
if(reyp > 1000.)
    cd0 = 0.44;
else
    cd0=24.0*(1.0+0.15*pow(reyp, 0.687))/reyp;
