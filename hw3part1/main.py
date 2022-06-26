from scipy.stats import norm
from matplotlib import pyplot as plt
import numpy as np
from vaccination import Vaccination #required ones are imported

my_vaccination=Vaccination() #equiated

input_Fuzzy_Set=[]#input fuzzy set's vector is defined as empty set
input_Fuzzy_Set.append([-0.00001,0,0.6])#first fuzzy set is appended to our vector as formulation [lower bound, mid bound, upper bound]
input_Fuzzy_Set.append([0.55,0.6,0.65])#second fuzzy set is appended
input_Fuzzy_Set.append([0.6,1,1.00001])#third fuzzy set is appended

output_Fuzzy_Set=[]#output fuzzy set's vector is defined as empty set
output_Fuzzy_Set.append([-0.25,-0.08,0])#first output fuzzy set is appended to our vector as formulation [lower bound, mid bound, upper bound]
output_Fuzzy_Set.append([-0.125,0,0.125])#second output fuzzy set is appended
output_Fuzzy_Set.append([0.025,0.125,0.25])#third output fuzzy set is appended

def calculate_membership(input,func,min,mid,max): #function to calculate membership is defined
    val=0#val is defined and set to zero
    if func==0:#since my memberships are in triangler form with these code its value is calculated as follows:
        if input<=mid:
            val= 1 / (mid-min) * (input-min)
        else:
            val= 1 / (mid - max) * (input -mid)+1
    else:
        val=norm.pdf(input, mid, max-min)/norm.pdf(mid,mid, max-min)
    if val>1:
        val=1
    elif val<0:
        val=0
    return  val

def calculate_set_membership(input,fuzzy_set,func):#function to set membership with set of value for each of set
    vals=[]# defined as empty array
    for set in fuzzy_set:
        vals.append(func(input,0,set[0],set[1],set[2]))#appended in required format
    return  vals

def calculate_mom(input_set,fuzzy_set,func):#function to calculate mean of maxima
    index_max=np.argmax(input_set) #highest value is calculated
    counter=0
    area=0
    i_begin=-1000
    i_ending=-1000

    for i in  np.linspace(fuzzy_set[index_max][0], fuzzy_set[index_max][2], 100):#for all i between max values for iterates
        val=func(i, 0, fuzzy_set[index_max][0], fuzzy_set[index_max][1], fuzzy_set[index_max][2])
        if val>=(input_set[index_max]*.99):
            counter=counter+1
            area=area+i
            if i_begin==-1000:
                i_begin=i
                i_ending=i
    return area/counter#gives mean of maxima value

def check_steady_state(m_array,error,set_point):
    counter=0
    index_max=np.argmax(m_array[0:int(len(m_array)/4)])
    print(index_max)
    for i in range(len(m_array)-1):
        diff=abs(m_array[i+1]-m_array[i])*50#difference is calculated and expanded with factor 50
        delta2=1
        if i>len(m_array)/2:
            delta2 = abs(np.average(m_array[i:-1] - set_point))*100
        if diff<error*0.025 and i>index_max and delta2<error:
            counter=counter+1
        if counter>5:
            return  i#after 5 case satisfies condition we return the 10 times day of reaching steady state (10x 3.4 =34 for example as i)
vaccination_array=[]
delta=0.25#define as max output value
for p in range(200):#for p from 0 to 200(for 20 days)#apply vaccination:
    vaccination,rate=my_vaccination.checkVaccinationStatus()
    inp_fuzzy=np.array(calculate_set_membership(vaccination, input_Fuzzy_Set, calculate_membership))
    inp_fuzzy=inp_fuzzy[::-1] # Adjust the fuzzy sets.
    delta = calculate_mom(inp_fuzzy, output_Fuzzy_Set, calculate_membership)

    vaccination_array.append(vaccination)
    print(p,vaccination,delta)
    my_vaccination.vaccinatePeople(delta)#applied

average=np.average(my_vaccination.vaccinated_percentage_curve_[100:200])#avarage is calculated
print("Overshoot:",np.max(my_vaccination.vaccinated_percentage_curve_)/average*100)#overshoot is displayed where displayed value=100+(%overshoot)
index=check_steady_state(my_vaccination.vaccinated_percentage_curve_,2,average)#ten times day of reaching steady state is found and equited
my_vaccination.viewVaccination(index,np.sum(my_vaccination.vaccination_rate_curve_[0:index]))#viewVaccination function of vaccination.py is applied to get visuals

fuzzy_set_array1=[]#all these sets are defined to store values for 10000 iteration
fuzzy_set_array2=[]
fuzzy_set_array3=[]
fuzzy_set_array4=[]
fuzzy_set_array5=[]
fuzzy_set_array6=[]

x = np.linspace(-1, 1, 10000)
for i in x:
    fuzzy_set_array1.append(calculate_membership(i, 0, input_Fuzzy_Set[0][0],
                                                 input_Fuzzy_Set[0][1], input_Fuzzy_Set[0][2]))#first three are for input values
    fuzzy_set_array2.append(calculate_membership(i, 0, input_Fuzzy_Set[1][0],
                                                 input_Fuzzy_Set[1][1], input_Fuzzy_Set[1][2]))
    fuzzy_set_array3.append(calculate_membership(i, 0, input_Fuzzy_Set[2][0],
                                                 input_Fuzzy_Set[2][1], input_Fuzzy_Set[2][2]))

    fuzzy_set_array4.append(calculate_membership(i, 0, output_Fuzzy_Set[0][0],
                                                 output_Fuzzy_Set[0][1], output_Fuzzy_Set[0][2]))#other three are for ourput values
    fuzzy_set_array5.append(calculate_membership(i, 0, output_Fuzzy_Set[1][0],
                                                 output_Fuzzy_Set[1][1], output_Fuzzy_Set[1][2]))
    fuzzy_set_array6.append(calculate_membership(i, 0, output_Fuzzy_Set[2][0],
                                                 output_Fuzzy_Set[2][1], output_Fuzzy_Set[2][2]))

test1=calculate_set_membership(0.6,input_Fuzzy_Set,calculate_membership)
test2=calculate_set_membership(0.51,input_Fuzzy_Set,calculate_membership)

plt.plot(x, fuzzy_set_array1)#plotting codes and arrangements
plt.plot(x, fuzzy_set_array2)
plt.plot(x, fuzzy_set_array3)
plt.grid()
plt.xticks(np.arange(-1.1, 1.1, 0.1))
plt.show()

plt.plot(x, fuzzy_set_array4)
plt.plot(x, fuzzy_set_array5)
plt.plot(x, fuzzy_set_array6)
plt.grid()
plt.xticks(np.arange(-1.1, 1.1, 0.1))
plt.show()

