int main(int argc, char *argv[])
{
    // MAIN LOOP
    for (int step = 2; step <= num_ts; step++)
    {

        if (debug)
        {
            printf("Step %d", step);
        }
        
        if (debug)
        {
            printf("Loop 1");
        }
        
        // Setup X, Y, Z, r and b from the updated initial conditions
        // LOOP 1
        for (int i = 1; i <= num_bodies; i++)
        {
            // LOOP 2
            for(int j = 1; j <= i-1; j++) // Check this
            {
                X[i][j][0] = x[j][0] - x[i][0];
                Y[i][j][0] = y[j][0] - y[i][0];
                Z[i][j][0] = z[j][0] - z[i][0];
                r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                b[i][j][0] = pow(r[i][j][0], -1.5);
            }

            // LOOP 3
            for (int j = i + 1; j < num_bodies; j++)
            {
                X[i][j][0] = x[j][0] - x[i][0];
                Y[i][j][0] = y[j][0] - y[i][0];
                Z[i][j][0] = z[j][0] - z[i][0];
                r[i][j][0] = pow(X[i][j][0], 2) + pow(Y[i][j][0], 2) + pow(Z[i][j][0], 2);
                b[i][j][0] = pow(r[i][j][0], -1.5);
            }
        }

        // LOOP 4
        for (int k = 1; k <= mac_degree; k++)
        {
            // LOOP 5
            for (int i = 1; i <= num_bodies; k++)
            {
                x[i][k] = u[i][k - 1]/k;
                y[i][k] = v[i][k - 1]/k;
                z[i][k] = w[i][k - 1]/k;
            }

            // LOOP 6
            for (int i = 1; i <= num_bodies; i++)
            {

                // LOOP 7
                for (int j = 0; j <= num_bodies; j++)
                {
                    X[i][j][k] = x[j][k] - x[i][k];
                    Y[i][j][k] = y[j][k] - y[i][k];
                    Z[i][j][k] = z[j][k] - z[i][k];
                    r[i][j][k] = cauchy_prod(X[i][j], X[i][j], k) + 
                                    cauchy_prod(Y[i][j], Y[i][j], k) + 
                                    cauchy_prod(Z[i][j], Z[i][j], k);
                    b[i][j][k] = cauchy_power(r[i][j], b[i][j], k - 1, -1.5);
                }
                
                // LOOP 8
                for (int j = i + 1; j <= num_bodies; j++)
                {
                    X[i][j][k] = x[j][k] - x[i][k];
                    Y[i][j][k] = y[j][k] - y[i][k];
                    Z[i][j][k] = z[j][k] - z[i][k];
                    r[i][j][k] = cauchy_prod(X[i][j], X[i][j], k) + 
                                    cauchy_prod(Y[i][j], Y[i][j], k) + 
                                    cauchy_prod(Z[i][j], Z[i][j], k);
                    b[i][j][k] = cauchy_power(r[i][j], b[i][j], k - 1, -1.5);
                }
            }
            
            // LOOP 9
            for (int i = 1; i <= num_bodies; i++)
            {
                u[i][k] = 0;
                v[i][k] = 0;
                w[i][k] = 0;
                
                // LOOP 10
                for (int j = 1; j <= i - 1; j++)
                {
                    u[i][k] += mass[j] * cauchy_prod(X[i][j], b[i][j], k - 1);
                    v[i][k] += mass[j] * cauchy_prod(Y[i][j], b[i][j], k - 1);
                    w[i][k] += mass[j] * cauchy_prod(Z[i][j], b[i][j], k - 1);
                }
                
                // LOOP 11
                for (int j = i + 1; j <= num_bodies; j++)
                {
                    u[i][k] += mass[j] * cauchy_prod(X[i][j], b[i][j], k - 1);
                    v[i][k] += mass[j] * cauchy_prod(Y[i][j], b[i][j], k - 1);
                    w[i][k] += mass[j] * cauchy_prod(Z[i][j], b[i][j], k - 1);
                }
                
                u[i][k] = u[i][k]/k;
                v[i][k] = v[i][k]/k;
                w[i][k] = w[i][k]/k;
            }
        }
        
        // Determine the values of the Maclaurin polynomial using Horner's algorithm and the
        // stored Maclauren coefficients
        for (int i = 1; i <= num_bodies; i++)
        {
            x_psm = horner_value(x[i], time_step, mac_degree);
            x[i][0] = x_psm;
            
            y_psm = horner_value(y[i], time_step, mac_degree);
            y[i][0] = y_psm;
            
            z_psm = horner_value(z[i], time_step, mac_degree);
            z[i][0] = z_psm;
            
            u_psm = horner_value(u[i], time_step, mac_degree);
            u[i][0] = u_psm;
            
            v_psm = horner_value(v[i], time_step, mac_degree);
            v[i][0] = v_psm;
            
            w_psm = horner_value(w[i], time_step, mac_degree);
            w[i][0] = w_psm;
        }

        // Output the step number based on the output granularity
        if (++num_cycles % granularity == 0)
        {
            printf("Step %d\n", step);
            
            if (verbose) 
            {
                for (int i = 1; i <= num_bodies; i++)
                {
                    printf("%lf\n",x[i][0]);
                    printf("%lf\n",y[i][0]);
                    printf("%lf\n",z[i][0]);
                }
            }
        } 
    }
}