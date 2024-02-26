#include <iostream>
#include <fftw3.h>

int main() {
    // Define the size of the input array
    const int N = 8;

    // Allocate memory for input and output arrays
    double in[N], out[N];

    // Initialize input array with sample data
    for (int i = 0; i < N; ++i) {
        in[i] = i + 1;
    }

    // Create FFTW plan for forward Fourier transform
    fftw_complex *fft_input = reinterpret_cast<fftw_complex*>(in);
    fftw_complex *fft_output = reinterpret_cast<fftw_complex*>(out);
    fftw_plan plan = fftw_plan_dft_1d(N, fft_input, fft_output, FFTW_FORWARD, FFTW_ESTIMATE);

    // Perform forward Fourier transform
    fftw_execute(plan);

    // Display results
    std::cout << "Input Data:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << in[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "FFT Output:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << out[i] << " ";
    }
    std::cout << std::endl;

    // Clean up FFTW resources
    fftw_destroy_plan(plan);
    fftw_cleanup();

    return 0;
}


//~ #include <iostream>
//~ #include <cmath>
//~ #include <fftw3.h>

//~ int main() {
    //~ const int N = 8; // Length of input signals
    //~ const int M = N * 2 - 1; // Length of output signal (convolution result)

    //~ // Define input signals
    //~ double input1[N] = {1, 2, 3, 4, 0, 0, 0, 0}; // Pad with zeros
    //~ double input2[N] = {1, 1, 1, 1, 0, 0, 0, 0}; // Pad with zeros
    //~ double result[M];

    //~ // Create FFTW plans
    //~ fftw_complex *fft_input1 = reinterpret_cast<fftw_complex*>(input1);
    //~ fftw_complex *fft_input2 = reinterpret_cast<fftw_complex*>(input2);
    //~ fftw_complex *fft_output = reinterpret_cast<fftw_complex*>(result);
    //~ fftw_plan forward_plan1 = fftw_plan_dft_1d(N, fft_input1, fft_input1, FFTW_FORWARD, FFTW_ESTIMATE);
    //~ fftw_plan forward_plan2 = fftw_plan_dft_1d(N, fft_input2, fft_input2, FFTW_FORWARD, FFTW_ESTIMATE);
    //~ fftw_plan inverse_plan = fftw_plan_dft_1d(M, fft_output, fft_output, FFTW_BACKWARD, FFTW_ESTIMATE);

    //~ // Perform forward FFTs
    //~ fftw_execute(forward_plan1);
    //~ fftw_execute(forward_plan2);

    //~ // Perform element-wise multiplication in frequency domain
    //~ for (int i = 0; i < N; ++i) {
        //~ fft_input1[i][0] *= fft_input2[i][0];
        //~ fft_input1[i][1] *= fft_input2[i][1];
    //~ }

    //~ // Perform inverse FFT to get convolution result
    //~ fftw_execute(inverse_plan);

    //~ // Normalize and print result
    //~ for (int i = 0; i < M; ++i) {
        //~ result[i] /= M;
        //~ std::cout << result[i] << " ";
    //~ }
    //~ std::cout << std::endl;

    //~ // Clean up FFTW resources
    //~ fftw_destroy_plan(forward_plan1);
    //~ fftw_destroy_plan(forward_plan2);
    //~ fftw_destroy_plan(inverse_plan);
    //~ fftw_cleanup();

    //~ return 0;
//~ }

