



def clamp(x, xmin, xmax):
    x = max(xmin, x)
    x = min(xmax, x)
    return x

def cap(ip, i, N):
    #return 3 + clamp(ip+i-3, -3, N+2) - ip 
    return 3 + min(ip+i-3, N+2) - ip 


print(-3, cap(-3,  0,  10))
print(-2, cap(-2,  0,  10))
print(-1, cap(-1,  0,  10))
print(0, cap( 0,  0,  10))
print(1, cap( 1,  0,  10))
print(2, cap( 2,  0,  10))
print(3, cap( 3,  0,  10))
print(4, cap( 4,  0,  10))
print(5, cap( 5,  6,  10))
print(6, cap( 6,  6,  10))
print(7, cap( 7,  6,  10))
print(8, cap( 8,  6,  10))
print(9, cap( 9,  6,  10))
print(10, cap( 10, 6,  10))
print(11, cap( 11, 6,  10))
print(12, cap( 12, 6,  10))
print(13, cap( 13, 6,  10))



