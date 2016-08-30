a=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
lb=4

for ub in range(15,25):
    lower=filter(lambda a:lb<a,a)
    mid=filter(lambda a:lb<a<=ub,a)
    upper=filter(lambda a:a>ub,a)

    print len(lower),len(mid),len(upper), len(lower)+len(mid)+len(upper)
