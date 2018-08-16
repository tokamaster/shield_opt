
import plotly
import plotly.offline as offline
from plotly.graph_objs import Scatter3d, Scatter, Layout

def rescale(OldValue,OldMin,OldMax,NewMin,NewMax):
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)
    NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    return NewValue

def next_prime():
    def is_prime(num):
        "Checks if num is a prime value"
        for i in range(2,int(num**0.5)+1):
            if(num % i)==0: return False
        return True

    prime = 3
    while(1):
        if is_prime(prime):
            yield prime
        prime += 2

def vdc(n, base=2):# van der Corput sequence
    vdc, denom = 0, 1
    while n:
        denom *= base
        n, remainder = divmod(n, base)
        vdc += remainder/float(denom)
    return vdc

def halton_sequence(size, dim):
    seq = []
    primeGen = next_prime()
    next(primeGen)
    for d in range(dim):
        base = next(primeGen)
        seq.append([vdc(i, base) for i in range(size)])
    return seq


number_of_datapoints=300
lower_x=-10
upper_x=10
lower_y=-10
upper_y=10
lower_z=-10
upper_z=10


points_to_search=[]
points_to_search_double_list = halton_sequence(number_of_datapoints, 3)
for x,y,z in zip(points_to_search_double_list[0],points_to_search_double_list[1],points_to_search_double_list[2]):
    new_x = rescale(x, 0.0, 1.0, lower_x, upper_x)
    new_y = rescale(y, 0.0, 1.0, lower_y, upper_y)
    new_z = rescale(z, 0.0, 1.0, lower_z, upper_z)



    points_to_search.append([new_x,new_y,new_z])
    #print(new_x,new_y,new_z)

#pip install plotly
#plotly.offline.plot({
#    "data": [Scatter3d(x=[item[0] for item in points_to_search], y=[item[1] for item in points_to_search],  z=[item[2] for item in points_to_search],mode='markers', marker=dict(size=3,color=list(range(0,len(points_to_search))) ,opacity=0.5,colorscale='Viridis'))],
#    "layout": Layout(title=str(len(points_to_search)) + " Halton sequence points for Steve")
#})
