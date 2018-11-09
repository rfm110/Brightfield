#!/usr/bin/env python
def generate_int(N):
    for i in xrange(N):
        yield i

def factorial(x):
    if x > 0:
        return x*factorial(x-1)
    else:
        return 1

def multiply2(generator):
    for element in generator:
        yield element*2
def sqrt(generator):
    for element in generator:
        yield element**2
a = generate_int(5)
# print len(list(a))
test = a
print repr(a)
o = multiply2(test)
q = sqrt(o)

# for i in o:
#     print i
# print q

for u in q:
    print u
