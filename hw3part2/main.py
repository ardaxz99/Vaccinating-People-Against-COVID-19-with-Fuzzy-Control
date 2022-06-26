from scipy.stats import norm
from matplotlib import pyplot as plt
import numpy as np
from vaccination import Vaccination#required ones are imported

my_vaccination=Vaccination()#equiated

input_Fuzzy_Set_current_rate=[]#input current states fuzzy set's vector is defined as empty set
input_Fuzzy_Set_current_rate.append([-0.00001,0,0.6])#first input current fuzzy set is appended to our vector as formulation [lower bound, mid bound, upper bound]
input_Fuzzy_Set_current_rate.append([0.58,0.6,0.62]))#second input current fuzzy set is appended
input_Fuzzy_Set_current_rate.append([0.6,1,1.00001]))#third input current fuzzy set is appended

input_Fuzzy_Set_failure_rate=[]#input failure states fuzzy set's vector is defined as empty set
input_Fuzzy_Set_failure_rate.append([-1 , -0.28 , 0.096])#first input failure fuzzy set is appended to our vector as formulation [lower bound, mid bound, upper bound]
input_Fuzzy_Set_failure_rate.append([-0.0921,0,0.0927])#second input failure fuzzy set is appended
input_Fuzzy_Set_failure_rate.append([0.000001,0.5,1.000001])#third input failure fuzzy set is appended

output_Fuzzy_Set=[]#output fuzzy set's vector is defined as empty set
output_Fuzzy_Set.append([-0.25  , -0.125,  -0.0625])#first output fuzzy set is appended to our vector as formulation [lower bound, mid bound, upper bound]
output_Fuzzy_Set.append([-0.125 ,- 0.0625,  -0.00001])#second output fuzzy set is appended
output_Fuzzy_Set.append([-0.0625,  0,  0.0625])#third output fuzzy set is appended
output_Fuzzy_Set.append([ 0.000001, 0.0625,  0.125])#fourth output fuzzy set is appended
output_Fuzzy_Set.append([0.0625,  0.125,  0.25])#fifth output fuzzy set is appended

rules=[]#rules vector is defined as empty set

rules.append(4) # 0
rules.append(3) # 1
rules.append(2) # 2
rules.append(1) # 3
rules.append(0) # 4

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

def calculate_set_membership(input_pi,input_pi_star,fuzzy_set1,fuzzy_set2,func):#function to set membership of two different input(input_Fuzzy_Set_failure_rate and input_Fuzzy_Set_current_rate) with set of values for each of set for each input
    vals1=[]# defined as empty array for input_Fuzzy_Set_current_rate
    vals2=[]# defined as empty array for input_Fuzzy_Set_failure_rate
    for set in fuzzy_set1:
        vals1.append(func(input_pi,0,set[0],set[1],set[2]))#appended in required format
    for set in fuzzy_set2:
        vals2.append(func(input_pi_star, 0, set[0], set[1], set[2]))#appended in required format
    return  vals1,vals2

def calculate_mom(input_set1,input_set2,fuzzy_set,_rules,func):#function to calculate mean of maxima

    index_max1=np.argmax(input_set1)#highest value is calculated for input_Fuzzy_Set_current_rate
    index_max2 = np.argmax(input_set2)#highest value is calculated for input_Fuzzy_Set_failure_rate
    rule_index=0
        #required rule to apply is found accordiing to which if will work:
    if index_max1==0 and index_max2==0: # Fuzzy Rule N
        rule_index=4
    elif index_max1==1 and index_max2==1: # Fuzzy Rule Z
        rule_index = 2
    elif index_max1==2 and index_max2==2: # Fuzzy Rule P
        rule_index = 0
    elif index_max1==0 and index_max2==1: # Fuzzy Rule SZ
        rule_index=3
    elif index_max1==0 and index_max2==2: # Fuzzy Rule SP
        rule_index=2
    elif index_max1==1 and index_max2==2: # Fuzzy Rule SP
        rule_index=1
    elif index_max1==1 and index_max2==0: # Fuzzy Rule SN
        rule_index=2
    elif index_max1==2 and index_max2==0: # Fuzzy Rule SP
        rule_index=1
    elif index_max1==2 and index_max2==1: # Fuzzy Rule P
        rule_index=0

    mu=(input_set1[index_max1]+input_set2[index_max2])/2

    if mu > 1:
        mu = 1
    counter = 0.00000001
    i_begin = -1000
    i_end = -1000
    print("R:", rule_index, index_max1, index_max2, input_set1, input_set2)
    area = 0
    for i in np.linspace(fuzzy_set[rule_index][0], fuzzy_set[rule_index][2], 100):  # Calculate the mom
        val = func(i, 0, fuzzy_set[rule_index][0], fuzzy_set[rule_index][1], fuzzy_set[rule_index][2])
        if val >= (mu * .95):
            counter = counter + 1
            area = area + i
            if i_begin == -1000:
                i_begin = i
            i_end = i

    return area / counter#returned as mean of maxima
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
        #if counter>5:
            return  i#i is returned as ten times day of reaching steady state (10x 3.4 =34 for example as i)
vaccination_array=[]
delta=0.25#define as max output value

for p in range(200):#for p from 0 to 200(for 20 days)#apply vaccination:
    vaccination,rate=my_vaccination.checkVaccinationStatus()
    inp_fuzzy1,inp_fuzzy2=np.array(calculate_set_membership(vaccination,rate, input_Fuzzy_Set_current_rate,
                                                            input_Fuzzy_Set_failure_rate, calculate_membership))

    delta = calculate_mom(inp_fuzzy1,inp_fuzzy2, output_Fuzzy_Set,rules, calculate_membership)

    vaccination_array.append(vaccination)
    print(p,vaccination,rate,delta)
    my_vaccination.vaccinatePeople(delta)#applied


average=np.average(my_vaccination.vaccinated_percentage_curve_[100:200])#avarage is calculated
print("Overshoot:",np.max(my_vaccination.vaccinated_percentage_curve_)/average*100)#overshoot is displayed where displayed value=100+(%overshoot)
index=check_steady_state(my_vaccination.vaccinated_percentage_curve_,2,average)#ten times day of reaching steady state is found and equited
if index<5:
    index=10

my_vaccination.viewVaccination(index,np.sum(my_vaccination.vaccination_rate_curve_[0:index]))#viewVaccination function of vaccination.py is applied to get visuals

fuzzy_set_array1=[]#all these sets are defined to store values for 10000 iteration for input_Fuzzy_Set_current_rate
fuzzy_set_array2=[]
fuzzy_set_array3=[]
fuzzy_set_array12=[]#all these sets are defined to store values for 10000 iteration for input_Fuzzy_Set_failure_rate
fuzzy_set_array22=[]
fuzzy_set_array32=[]

fuzzy_set_array4=[]#all these sets are defined to store values for 10000 iteration for output_Fuzzy_Set
fuzzy_set_array5=[]
fuzzy_set_array6=[]
fuzzy_set_array7=[]
fuzzy_set_array8=[]



x = np.linspace(-1, 1, 10000)
for i in x:
    fuzzy_set_array1.append(calculate_membership(i, 0, input_Fuzzy_Set_current_rate[0][0],
                                                 input_Fuzzy_Set_current_rate[0][1], input_Fuzzy_Set_current_rate[0][2]))#first three are for input_Fuzzy_Set_current_rate values
    fuzzy_set_array2.append(calculate_membership(i, 0, input_Fuzzy_Set_current_rate[1][0],
                                                 input_Fuzzy_Set_current_rate[1][1], input_Fuzzy_Set_current_rate[1][2]))
    fuzzy_set_array3.append(calculate_membership(i, 0, input_Fuzzy_Set_current_rate[2][0],
                                                 input_Fuzzy_Set_current_rate[2][1], input_Fuzzy_Set_current_rate[2][2]))

    fuzzy_set_array12.append(calculate_membership(i, 0, input_Fuzzy_Set_failure_rate[0][0], input_Fuzzy_Set_failure_rate[0][1],#These three are for input_Fuzzy_Set_failure_rate values
                                      input_Fuzzy_Set_failure_rate[0][2]))
    fuzzy_set_array22.append(calculate_membership(i, 0, input_Fuzzy_Set_failure_rate[1][0], input_Fuzzy_Set_failure_rate[1][1],
                                      input_Fuzzy_Set_failure_rate[1][2]))
    fuzzy_set_array32.append(calculate_membership(i, 0, input_Fuzzy_Set_failure_rate[2][0], input_Fuzzy_Set_failure_rate[2][1],
                                      input_Fuzzy_Set_failure_rate[2][2]))

    fuzzy_set_array4.append(calculate_membership(i, 0, output_Fuzzy_Set[0][0],
                                                 output_Fuzzy_Set[0][1], output_Fuzzy_Set[0][2]))#These three are for output_Fuzzy_Set values
    fuzzy_set_array5.append(calculate_membership(i, 0, output_Fuzzy_Set[1][0],
                                                 output_Fuzzy_Set[1][1], output_Fuzzy_Set[1][2]))
    fuzzy_set_array6.append(calculate_membership(i, 0, output_Fuzzy_Set[2][0],
                                                 output_Fuzzy_Set[2][1], output_Fuzzy_Set[2][2]))
    fuzzy_set_array7.append(calculate_membership(i, 0, output_Fuzzy_Set[3][0],
                                                 output_Fuzzy_Set[3][1], output_Fuzzy_Set[3][2]))
    fuzzy_set_array8.append(calculate_membership(i, 0, output_Fuzzy_Set[4][0],
                                                 output_Fuzzy_Set[4][1], output_Fuzzy_Set[4][2]))

plt.plot(x, fuzzy_set_array1)#plotting codes and arrangements
plt.plot(x, fuzzy_set_array2)
plt.plot(x, fuzzy_set_array3)
plt.grid()
plt.xticks(np.arange(-1.1, 1.1, 0.1))
plt.show()

plt.plot(x, fuzzy_set_array12)
plt.plot(x, fuzzy_set_array22)
plt.plot(x, fuzzy_set_array32)
plt.grid()
plt.xticks(np.arange(-1.1, 1.1, 0.1))
plt.show()

plt.plot(x, fuzzy_set_array4)
plt.plot(x, fuzzy_set_array5)
plt.plot(x, fuzzy_set_array6)
plt.plot(x, fuzzy_set_array7)
plt.plot(x, fuzzy_set_array8)
plt.grid()
plt.xticks(np.arange(-1.1, 1.1, 0.1))
plt.show()

