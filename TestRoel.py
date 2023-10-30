import scipy.sparse as sp

# Define the dimensions of the identity matrix
Nx = 5  # Change this to your desired size

# Create the sparse identity matrix in CSR format
PHI = sp.identity(Nx, dtype='float', format='csr')

# Print the sparse identity matrix
print(PHI)
