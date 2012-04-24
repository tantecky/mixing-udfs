Thread* thread_l;
Thread* thread_s;

real slip_x;
real slip_y;
real slip_z;
real slip;

real reyp;
real cd0;

/*liquid phase - primary*/
thread_l = THREAD_SUB_THREAD(mix_thread, s_col);
/*solid phase - secondary*/
thread_s = THREAD_SUB_THREAD(mix_thread, f_col);

slip_x = C_U(cell, thread_l) - C_U(cell, thread_s);
slip_y = C_V(cell, thread_l) - C_V(cell, thread_s);
slip_z = C_W(cell, thread_l) - C_W(cell, thread_s);

/*Euclidean norm of slip velocity*/
slip = sqrt(slip_x*slip_x + slip_y*slip_y + slip_z*slip_z);

/*relative Reynolds number for the primary phase*/
reyp = RHO_L*slip*DIAMETER/MU_L;

/*drag coefficient*/
/*Schiller-Nauman*/
if(reyp > 1000.)
    cd0 = 0.44;
else
    cd0=24.0*(1.0+0.15*pow(reyp, 0.687))/reyp;
