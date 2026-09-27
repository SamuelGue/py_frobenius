import math
import numpy as np

class Frobenius:
    def __init__(self, *args):
        self.args = sorted(args)

        def validate_args(self):
            if len(self.args) < 2:
                raise ValueError("At least two arguments are required.")
            elif any(not isinstance(arg, int) or arg <= 0 for arg in self.args):
                raise ValueError("All arguments must be positive integers.")
            elif math.gcd(*self.args) != 1:
                raise ValueError("Arguments must be coprime.")
            else:
                return True

        validate_args(self)

    def frobenius_number(self):
        # If there are only two arguments we use Sylvester's formula for the Frobenius number
        # If there are more than two arguments we use the BFD algorithm
        if len(self.args) == 2:
            a, b = self.args
            return a * b - a - b
        else:
            A = np.array(self.args).astype(int)
            Q = [0]
            P = np.zeros(A[0]).astype(int)
            P[0] = len(self.args)
            S = np.ones(A[0]).astype(int) * A[0] * A[-1]
            S[0] = 0
            Amod = A % A[0]

            while Q:
                v = Q.pop(0)
                for j in range(1, P[v]):
                    u = v + Amod[j]
                    if u >= A[0]:
                        u -= A[0]
                    u = int(u)
                    w  = S[v] + A[j]
                    w = int(w)
                    if w < S[u]:
                        S[u] = w
                        P[u] = j+1
                        if u not in Q:
                            Q.append(u)

            return int(max(S) - A[0])

    def find_representations(self, n, out_dict=None):

        if out_dict is None:
            out_dict = {i: 0 for i in self.args}
        
        if n in self.args:
            result = dict(out_dict)
            result[n] += 1
            return [result]

        if n < self.args[0]:
            return None

        results = []
        for arg in self.args:
            sub = self.find_representations(n-arg, {**out_dict, arg: out_dict[arg]+1})
            if sub is not None:
                results.extend(sub)

        return [dict(t) for t in {tuple(d.items()) for d in results}] if results else None

    def denumerant(self, N, up_to_N = False):
        # We use the Cooley-Tukey Fast Fourier Transform (FFT) algorithm to compute the product of polynomials. To compute the denumerant of N, we determine the Nth coefficient of the generating function of the denumerants.
        def recursive_FFT(a):
            # Uses the Fast Fourier Transform to evaluate a polynomial 'a' at the roots of unity
            n = len(a)

            if n == 1:
                return [a[0]]

            theta = -2 * math.pi / n
            w = list(complex(math.cos(theta * i), math.sin(theta * i)) for i in range(n))

            a_E = a[0::2]
            a_O = a[1::2]

            y_E = recursive_FFT(a_E)
            y_O = recursive_FFT(a_O)

            y = [0] * n

            for k in range(0, n//2):
                y[k] = y_E[k] + w[k]*y_O[k]
                y[k+n//2] = y_E[k] - w[k]*y_O[k]

            return y

        def inverse_recursive_FFT(a):
            # Computes the coefficients of a polynomial with values 'a' at the roots of unity
            n = len(a)

            if n == 1:
                return [a[0]]

            theta = -2 * math.pi / n
            w = list(complex(math.cos(theta * i), -math.sin(theta * i)) for i in range(n))

            a_E = a[0::2]
            a_O = a[1::2]

            y_E = inverse_recursive_FFT(a_E)
            y_O = inverse_recursive_FFT(a_O)

            y = [0] * n

            for k in range(0, n//2):
                y[k] = y_E[k] + w[k]*y_O[k]
                y[k+n//2] = y_E[k] - w[k]*y_O[k]

            return y

        def poly_prod(f,g):
            # Computes the product of polynomials f and g. Cooley-Tukey requires polynomials wiht degree equal to a power of 2, so we pad the polynomials with zeros.
            result_len = len(f) + len(g) - 1

            N = 1
            while N < result_len:
                N *= 2

            f_padded = f + [0] * (N - len(f))
            g_padded = g + [0] * (N - len(g))

            a = recursive_FFT(f_padded)
            b = recursive_FFT(g_padded)

            c = [a[i] * b[i] for i in range(N)]

            y = inverse_recursive_FFT(c)
            result = [round(i.real / N) for i in y]

            return result[:result_len]

        def make_poly(arg, n=N):
            # Computes the required terms of the expression 1/(1-x^a_i)
            highest_power = n // arg
            coeffs = np.zeros(1 + highest_power * arg)
            for i in range(0,highest_power + 1):
                coeffs[i * arg] = 1
            return coeffs.astype(int).tolist()

        poly = make_poly(self.args[0],n=N)
        for i in range(1,len(self.args)):
            poly = poly_prod(poly,make_poly(self.args[i],n=N))

        if up_to_N:
            return poly[:N]
        else:
            return poly[N]

    def unrepresented(self):
        # Returns all integers that have no representation
        frob = self.frobenius_number()
        unreps = []
        denumerant_func = self.denumerant(frob+1, up_to_N=True)
        unreps = [i for i in range(0,frob + 1) if denumerant_func[i] == 0]
        return unreps

    def num_unrepresented(self):
        return len(self.unrepresented())
