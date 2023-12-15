# How to get the final $y$ from the experiment interferogram

I believe I have a solution on how to get the interferogram associated to the
sparse spectra from the experimental interfegram acquaried the Inbuia line. 

From our talk today what they do today is:

1. They acquire an interferogram with the background information, let us call it
   $I_B$. This can be done "once".

1. They then acquire interferograms from experimental samples, let us denote it
   $I_E$.

1. Then compute the Fourier transform of boh interferogram moving to the
   frequency domain using the Fourier transform. That is they compute $y_B =
   \mathcal{F}[I_B]$, and $y_E = \mathcal{F}[I_E]$.

1. The final interfegram, the one that is supposed to be sparse is $y = y_E /
   y_B$.

Great, now I will show conceptually (I cannot do the details) how to get the
final $y$ as a Fourier transform of the original interfegrams. Start by
computing 
$$
    I_B^{-1} = \mathcal{F}^{-1}[1 / y_B].
$$

From the formulas above we have
$$
y = \mathcal{F}[I_E] \cdot \mathcal{F}[I_B^{-1}].
$$
Now we can use the convolution theorem and conclude that
$$
\mathcal{F}^{-1}[y] = I_E \star I_B^{-1}.
$$
Done.
