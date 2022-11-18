n = 100


M = []


for i in range(n * n):
    M.append(0)

for i in range(n):
    M[n * i + i] = 2
    for j in range(5):
        M[n * i + ((i * + j * (n - 13)) % n)] = 1
        M[n * i + ((i * - j * (n - 7)) % n)] = 1

b = []

for i in range(n):
    b.append(i)
print("M:")
for i in range(n):
    for j in range(n):
        print(M[i * n + j], end=", ")
    print()
print()
print("b:")
for x in b:
    print(x, end=", ")
print()

Mb = []

for i in range(n):
    elem = 0
    for j in range(n):
        elem += M[i * n + j] * b[j]
    Mb.append(elem)
print("Mb:")
for x in Mb:
    print(x, end=", ")
print()

sparse_M = [[] for _ in range(n)]
indecie_rows = [[] for _ in range(n)]
for i in range(n):
    for j in range(n):
        if M[i * n + j] != 0:
            indecie_rows[i].append(j)
            sparse_M[i].append(M[i * n + j])


max_number_of_nonzero_elements_on_row = 32

print("Sparse M:")
for row in sparse_M:
    if len(row) < max_number_of_nonzero_elements_on_row:
        for x in row:
            print(x, end=", ")
        for _ in range(max_number_of_nonzero_elements_on_row - len(row)):
            print("0", end=", ")
    else:
        print("ERROR! More elements on row than there should")
    print()
print()

print("Indecies of sparse M:")
for row in indecie_rows:
    if len(row) < max_number_of_nonzero_elements_on_row:
        for x in row:
            print(x, end=", ")
        for _ in range(max_number_of_nonzero_elements_on_row - len(row)):
            print("0", end=", ")
    else:
        print("ERROR! More elements on row than there should")
    print()
print()
