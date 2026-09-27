from py_frobenius import Frobenius

frob = Frobenius(3,5,7)

print(f"Frobenius Number: {frob.frobenius_number()}")

print(f"Representations of 21: {frob.find_representations(21)}")

print(f"Number of representations of 21: {frob.denumerant(21)}")

print(f"All unrepresented integers: {frob.unrepresented()}")

print(f"Number of unrepresented integers: {frob.num_unrepresented()}")