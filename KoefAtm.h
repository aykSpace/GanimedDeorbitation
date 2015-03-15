double C[36][10] = {
        {-0.042e0, 0.001e0, 0.024e0, 0.068e0,-0.006e0, 0.002e0, 0.123e0, 0.059e0, 0.014e0,-0.014e0},
        {-0.04e0 , 0.e0   , 0.e0   , 0.023e0, 0.002e0,-0.003e0, 0.031e0, 0.022e0, 0.003e0,-0.008e0},
        {-0.032e0,-0.002e0, 0.e0   , 0.014e0, 0.002e0,-0.002e0, 0.001e0, 0.009e0, 0.001e0,-0.003e0},
        {-0.019e0,-0.004e0, 0.006e0, 0.006e0, 0.002e0,-0.002e0,-0.019e0,-0.003e0,-0.002e0, 0.002e0},
        { 0.008e0, 0.001e0, 0.009e0,-0.005e0,-0.001e0,-0.006e0,-0.072e0,-0.014e0,-0.003e0, 0.01e0 },
        { 0.064e0, 0.002e0,-0.005e0,-0.027e0,-0.002e0,-0.005e0,-0.104e0,-0.015e0,-0.003e0, 0.012e0},
        { 0.169e0, 0.001e0, 0.005e0,-0.064e0,-0.002e0, 0.009e0,-0.086e0,-0.039e0,-0.002e0, 0.012e0},
        { 0.093e0, 0.001e0, 0.003e0,-0.066e0, 0.002e0, 0.008e0,-0.069e0,-0.07e0, -0.002e0, 0.008e0},
        { 0.031e0, 0.002e0, 0.006e0,-0.066e0, 0.004e0, 0.007e0,-0.05e0 ,-0.116e0, 0.004e0, 0.002e0},
        {-0.001e0, 0.003e0,-0.001e0,-0.071e0, 0.004e0, 0.009e0,-0.046e0,-0.154e0, 0.025e0,-0.001e0},
        {-0.009e0, 0.003e0,-0.004e0,-0.076e0, 0.003e0, 0.01e0 ,-0.049e0,-0.172e0, 0.036e0,-0.002e0},
        {-0.003e0, 0.005e0,-0.012e0,-0.08e0 ,-0.009e0, 0.017e0,-0.067e0,-0.211e0, 0.06e0 ,-0.002e0},
        { 0.019e0, 0.005e0,-0.004e0,-0.094e0,-0.001e0, 0.018e0,-0.07e0, -0.241e0, 0.078e0, 0.001e0},
        { 0.036e0, 0.002e0, 0.002e0,-0.108e0, 0.007e0, 0.018e0,-0.073e0,-0.267e0, 0.095e0, 0.005e0},
        { 0.051e0,-0.002e0, 0.006e0,-0.12e0 , 0.014e0, 0.019e0,-0.078e0,-0.293e0, 0.109e0, 0.008e0},
        { 0.068e0,-0.013e0, 0.009e0,-0.128e0, 0.02e0 , 0.021e0,-0.083e0,-0.316e0, 0.122e0, 0.012e0},
        { 0.082e0,-0.019e0, 0.009e0,-0.13e0 , 0.026e0, 0.027e0,-0.094e0,-0.337e0, 0.133e0, 0.015e0},
        { 0.086e0,-0.02e0 , 0.007e0,-0.134e0, 0.03e0 , 0.026e0,-0.099e0,-0.356e0, 0.141e0, 0.019e0},
        { 0.084e0,-0.019e0, 0.008e0,-0.138e0, 0.037e0, 0.025e0,-0.118e0,-0.391e0, 0.149e0, 0.031e0},
        { 0.081e0,-0.019e0, 0.005e0,-0.152e0, 0.043e0, 0.024e0,-0.128e0,-0.419e0, 0.171e0, 0.033e0},
        {0.087e0 ,-0.018e0, 0.001e0,-0.165e0, 0.048e0, 0.024e0,-0.125e0,-0.433e0, 0.178e0, 0.035e0},
        {0.097e0 ,-0.022e0,-0.003e0,-0.175e0, 0.052e0, 0.024e0,-0.119e0,-0.444e0, 0.184e0, 0.035e0},
        {0.102e0 ,-0.023e0,-0.013e0,-0.184e0, 0.053e0, 0.029e0,-0.113e0,-0.452e0, 0.187e0, 0.034e0},
        {0.105e0 ,-0.023e0,-0.02e0 ,-0.186e0, 0.054e0, 0.034e0,-0.102e0,-0.453e0, 0.188e0, 0.034e0},
        {0.107e0 ,-0.02e0 ,-0.027e0,-0.189e0, 0.054e0, 0.033e0,-0.088e0,-0.448e0, 0.185e0, 0.033e0},
        {0.102e0 ,-0.011e0,-0.02e0 ,-0.181e0, 0.053e0, 0.032e0,-0.066e0,-0.436e0, 0.172e0, 0.034e0},
        {0.059e0 ,-0.002e0,-0.016e0,-0.016e0, 0.051e0, 0.028e0,-0.045e0,-0.42e0 , 0.158e0, 0.031e0},
        {0.019e0 , 0.007e0,-0.006e0,-0.125e0, 0.049e0, 0.029e0,-0.023e0,-0.393e0, 0.145e0, 0.026e0},
        {-0.023e0, 0.015e0,-0.004e0,-0.077e0, 0.045e0, 0.015e0, 0.008e0,-0.344e0, 0.129e0, 0.013e0},
        {0.008e0 , 0.021e0, 0.015e0, 0.008e0, 0.036e0, 0.017e0, 0.039e0,-0.24e0 , 0.11e0 ,-0.02e0 },
        {0.051e0 ,-0.021e0, 0.01e0 , 0.123e0, 0.026e0,-0.023e0, 0.073e0,-0.105e0, 0.091e0,-0.069e0},
        {0.075e0 ,-0.057e0, 0.033e0, 0.214e0, 0.018e0,-0.022e0, 0.153e0,-0.034e0, 0.075e0,-0.094e0},
        {0.105e0 ,-0.07e0 , 0.074e0, 0.24e0 , 0.007e0,-0.021e0, 0.229e0, 0.054e0, 0.053e0,-0.101e0},
        {0.137e0 ,-0.072e0, 0.131e0, 0.243e0,-0.002e0,-0.027e0, 0.272e0, 0.056e0, 0.017e0,-0.07e0 },
        {0.179e0 ,-0.079e0, 0.166e0, 0.207e0,-0.008e0,-0.035e0, 0.265e0, 0.024e0, 0.023e0,-0.046e0},
        {0.228e0 ,-0.09e0 , 0.19e0 , 0.174e0,-0.012e0,-0.041e0, 0.245e0,-0.014e0, 0.024e0,-0.03e0 },
};

double F[29][10] = {
        {  1.5e0  ,- 0.5e0 ,-0.4e0 ,-0.1e0 ,- 0.8e0,- 0.3e0,  4.3e0,  0.5e0,  2.e0 ,  0.3e0},
        {- 0.8e0  ,- 3.e0  , 2.8e0 , 0.5e0 ,- 1.2e0,  0.2e0,  6.3e0, -0.6e0,  1.3e0,  0.3e0},
        {  0.e0   ,- 5.6e0 , 5.8e0 , 0.4e0 ,  0.1e0,- 0.2e0,  8.3e0, -0.8e0,  0.6e0,  0.5e0},
        {  0.8e0  ,- 8.1e0 , 7.9e0 , 0.1e0 ,  0.1e0,  0.e0 , 10.4e0, -1.e0 ,- 0.2e0,  0.6e0},
        {  1.e0   ,- 9.e0  , 15.8e0,-4.e0  ,  1.2e0,  3.2e0,  9.e0 , -1.7e0,  0.1e0,  0.1e0},
        {  0.4e0  ,- 5.4e0 , 12.2e0, 0.5e0 ,  0.e0 ,  0.2e0,  6.1e0, -1.1e0,  1.3e0,- 0.1e0},
        {- 1.6e0  ,- 3.3e0 , 11.e0 , 4.6e0 ,- 1.6e0,- 1.4e0,  5.2e0, -0.2e0,  1.2e0,- 1.e0 },
        {- 4.3e0  ,- 1.8e0 , 10.7e0, 9.8e0 ,- 2.3e0,- 2.1e0,  5.1e0,  1.3e0,  1.3e0,- 1.7e0},
        {- 7.8e0  ,- 1.1e0 , 11.2e0, 13.4e0,- 3.3e0,- 2.6e0,  5.4e0,  3.3e0,  0.8e0,- 2.4e0},
        {-11.4e0  ,- 0.6e0 , 12.4e0, 16.6e0,- 4.3e0,- 4.2e0,  5.9e0,  5.4e0,  0.4e0,- 3.1e0},
        {-20.e0   ,  0.e0  , 13.9e0, 19.6e0,- 4.8e0,- 6.4e0,  6.3e0,  7.5e0,  0.e0 ,- 3.8e0},
        {-17.5e0  ,- 2.5e0 , 15.8e0, 26.9e0,- 6.e0 ,- 6.9e0,  8.5e0, 13.3e0,- 0.4e0,- 4.4e0},
        {-15.e0   ,-15.e0  , 17.9e0, 33.8e0,- 7.4e0,- 7.2e0, 10.8e0, 19.1e0,- 0.8e0,- 5.e0 },
        {-12.5e0  ,-27.5e0 , 19.8e0, 41.2e0,- 8.6e0,- 7.7e0, 13.1e0, 25.e0 ,- 1.3e0,- 5.6e0},
        {- 5.e0   ,-35.e0  , 22.e0 , 48.9e0,-10.e0 ,- 8.e0 , 15.8e0, 30.8e0,- 2.5e0,- 6.7e0},
        {  5.9e0  ,-24.1e0 , 20.5e0, 51.e0 ,-10.e0 ,- 9.5e0, 16.7e0, 32.9e0,- 6.3e0,- 9.6e0},
        { 16.7e0  ,-13.3e0 , 19.2e0, 53.5e0,-10.e0 ,-10.8e0, 17.5e0, 35.e0 ,-10.e0 ,-12.5e0},
        { 27.5e0  ,- 2.5e0 , 16.4e0, 53.5e0,- 7.4e0,- 7.9e0, 15.4e0, 32.5e0,-11.7e0,-12.9e0},
        { 17.5e0  ,  7.5e0 , 10.8e0, 48.5e0,- 5.2e0,- 7.4e0, 13.3e0, 30.e0, -13.3e0,-13.3e0},
        {  7.5e0  ,17.5e0  ,  9.2e0, 36.e0 ,- 2.5e0,- 3.2e0, 10.5e0, 24.e0, -10.e0 ,- 9.5e0},
        {- 1.6e0  , 8.3e0  , 12.5e0, 27.e0 ,- 5.e0 ,  2.5e0, 12.3e0, 20.5e0,-10.e0 ,- 7.8e0},
        {-10.8e0  ,-0.8e0  , 15.5e0, 24.e0 ,- 4.e0 ,  8.5e0, 14.e0 , 17.e0 ,-10.e0 ,- 6.e0 },
        {-10.e0   ,-20.e0  , 19.3e0, 12.1e0,-12.5e0, 16.8e0, 18.8e0, 16.7e0,-10.e0 ,- 4.1e0},
        {  0.e0   ,-30.e0  , 25.8e0,  7.1e0,-14.6e0, 21.2e0, 23.6e0, 16.4e0,-10.e0 ,- 2.1e0},
        {  6.e0   ,-24.e0  , 30.7e0, 10.2e0,-12.2e0, 19.5e0, 28.4e0, 16.1e0,-10.e0 ,- 0.2e0},
        {  5.8e0  ,-11.8e0 , 32.8e0, 19.8e0,- 8.3e0, 14.e0 , 33.1e0, 15.8e0,-10.e0 ,  1.7e0},
        {  5.5e0  ,  0.5e0 , 32.4e0, 28.9e0,-10.9e0, 11.5e0, 37.9e0, 15.5e0,-10.e0 ,  3.6e0},
        {  5.3e0  , 12.8e0 , 29.2e0, 38.4e0,-18.e0 , 11.2e0, 42.7e0, 15.3e0,-10.e0 ,  5.6e0},
        {  9.5e0  , 25.e0  , 26.2e0, 47.5e0,-25.e0 , 11.2e0, 47.5e0, 15.e0 ,-10.e0 ,  7.5e0},
};