
# coding: utf-8

# In[1]:


def extract_last2bits(n):
    x=n>>2
    y=x<<2
    return(n-y)


# In[2]:


def extract_lastkbits(n,k):
    x=n>>k
    y=x<<k
    return(n-y)


# In[3]:


#-------------------------------------------------------------------------
#Finds the next highest integer with a fixed number of 0's and 1's
#-------------------------------------------------------------------------


def next_highest_int(n): 
    c0=0
    c1=0
    c = n
    #Calculating c0
    while(((c&1)==0) and (c!=0)):
        c0=c0+1
        c=c>>1
    #Calculating c1
    while((c&1)==1):
        c1=c1+1
        c=c>>1
    #If there is no bigger number with the same no. of 1's
    if (c0 +c1 == 31 or c0 +c1== 0):
        return -1;
    #Calculating the position of righmost non-trailing 0
    p = c0 + c1;
    n |= (1 << p);
    #Clear all bits to the right of p
    n=n&~((1 << p) - 1);
    #Insert (c1-1) ones on the right.
    n=n|(1 << (c1 - 1)) - 1;
    return n;


# In[4]:


#--------------------------------------------------------------------------------------------------------------------------
#This function retuns the integers for the states in ascending order as well as their occupancies in two separate arrays.
#n is the number of orbitals and e is the number of electrons
#--------------------------------------------------------------------------------------------------------------------------


def states_with_MS0(n,e): 
    lowest=2**e-1
    z=lowest
    highest=(2**(2*n-e))*lowest
    MS0_states=[]
    all_orb_occu=[]
    while(z<=highest):
        x=z
        orb_occu=[]    #stores 0,1,2,3 depending on whether orbital has 0, up, down or both electrons
        mstot=0
        
        #---------------------------------------------------------
        # Calculating MStotal
        #---------------------------------------------------------
        
        while(x>0):
            a=extract_last2bits(x)
            orb_occu.append(a)
            if(a==1):
                mstot=mstot+1
            if(a==2):
                mstot=mstot-1
            x=x>>2
        if(mstot==0):
            MS0_states.append(z)
            all_orb_occu.append(orb_occu) #stores orbital occupancy corresponding to a state.
        z=next_highest_int(z)
    return(MS0_states, all_orb_occu)


# In[5]:


#----------------------------------------------------------------------------------------------------------------------
# This function projects the states into either the ionic or covalent subspace
#
#----------------------------------------------------------------------------------------------------------------------

def ionic_covalent_space(n, e):
    all_states = states_with_MS0(n,e)[0]
    ionic = []
    covalent = []
    L = len(all_states)
    i = 0
    while(i<L):
        E = all_states[i]
        J = shiba(all_states[i], n)
        if(abs(E) != abs(J)):
            covalent.append([E,J])
            ionic.append([E,-J])
            all_states.remove(abs(J))
            L = len(all_states)
        if(E == J):
            covalent.append([E])
        if(E == -J):
            ionic.append([E])
        i = i + 1
    
    return(all_states, covalent, ionic)


# In[6]:


#----------------------------------------------------------------------------------------------------------------------
#Extracts the number of ones before the kth digit in z and return 1 or -1 depending on whether count is even or odd
#----------------------------------------------------------------------------------------------------------------------

def extract_phase(z,k):
    x=abs(z)
    count=0
    i=0
    while(i<k-1):
        count=count+(x&1)
        x=x>>1
        i=i+1
    if(count%2==0):
        return(1)
    else:
        return(-1)
 #   return(count)


# In[7]:


#------------------------------------------------------------------------------------------------------------------
# Creates up or down spin electron in the integer representation z at the kth orbital with the appropriate phase
# '10' for up electrons and '01' for down electrons
#------------------------------------------------------------------------------------------------------------------

def create_electron(z,k,spin): #Creates up or down spin electron in the integer representation z at the kth orbital
                               #'10' for up electrons and '01' for down electrons
    if (z==None):
        return None
    else:
        x=abs(z)
        initial_phase=int(z/x)
        a=x>>(2*(k-1))
        b=extract_last2bits(a)
        if(spin=='up'):
            phase=initial_phase*extract_phase(z,2*k)
            if((b==2)|(b==3)):
                return(None)
            else:
                return(phase*(x+(2**(2*k-1))))
        else:
            phase=initial_phase*extract_phase(z,2*k-1)
            if((b==1)|(b==3)):
                return(None)
            else:
                return(phase*(x+(2**(2*(k-1)))))


# In[8]:


#------------------------------------------------------------------------------------------------------------------
# Destroys up or down spin electron in the integer representation z at the kth orbital with the appropriate phase
# '10' for up electrons and '01' for down electrons
#------------------------------------------------------------------------------------------------------------------

def destroy_electron(z,k,spin):
    if (z==None):
        return None
    else:
        x=abs(z)
        initial_phase=int(z/x)
        a=x>>(2*(k-1))
        b=extract_last2bits(a)
        if(spin=='up'):
            phase=initial_phase*extract_phase(z,2*k)
            if((b==0)|(b==1)):
                return(None)
            else:
                return(phase*(x-(2**(2*k-1))))
        else:
            phase=initial_phase*extract_phase(z,2*k-1)
            if((b==0)|(b==2)):
                return(None)
            else:
                return(phase*(x-(2**(2*(k-1)))))


# In[9]:


#-------------------------------------------------------------------------------------------------------------------
# This function gives out the final states for a given state after hamiltonian acts on it as an array P
# Array A gives the absolute integer values of the final states
# Array count gives the sum of coefficients of the elements in array P
# X=integer value of the MS0 state, n=number of orbitals
#-------------------------------------------------------------------------------------------------------------------

import numpy as np

def final_states(x,n,t,U):
    i=1
    P=[]
    P_coeff=[]
    #--------------------------------------
    # Creating the array P
    #--------------------------------------
    
    while(i<n):
        a1=destroy_electron(x,i+1,'up')
        a2=create_electron(a1,i,'up')
        if(a2!=None):
            P.append(a2)
            P_coeff.append(t)
        b1=destroy_electron(x,i+1,'down')
        b2=create_electron(b1,i,'down')
        if(b2!=None):
            P.append(b2)
            P_coeff.append(t)
        c1=destroy_electron(x,i,'down')
        c2=create_electron(c1,i+1,'down')
        if(c2!=None):
            P.append(c2)
            P_coeff.append(t)
        d1=destroy_electron(x,i,'up')
        d2=create_electron(d1,i+1,'up')
        if(d2!=None):
            P.append(d2)
            P_coeff.append(t)
        i=i+1
        
    i=1
    while(i<=n):
        e1=create_electron(x,i,'up')
        e2=create_electron(e1,i,'down')
        e3=destroy_electron(e2,i,'down')
        e4=destroy_electron(e3,i,'up')
        if(e4!=None):
            P.append(e4)
            P_coeff.append(U)
        i=i+1
    
    
    b=0
    A=[]
    count=[]
    while(b<len(P)):
        j=0
        flag=0 
        #-------------------------------------------------------------------------
        # Flag variable to signal if an element has been encountered before
        # Checking if an element has occured before and updating count
        #-------------------------------------------------------------------------
        
        while(j<len(A)):
            if(abs(P[b])==A[j]):
                count[j]=count[j]+int(P[b]/abs(P[b]))*P_coeff[b]
                flag=1
            j=j+1
        
        #-------------------------------------------------------------------------
        # If the element has not been encountered, append it to the array A
        #-------------------------------------------------------------------------
        
        if(flag==0):
            A.append(abs(P[b]))
            count.append(int(P[b]/abs(P[b]))*P_coeff[b])
        b=b+1
    
#    P = np.asarray(P, dtype = int)
#    A = np.asarray(A, dtype = int)
#    count = np.asarray(count, dtype = np.float64)
        
    return(P,A,count)


# In[10]:


#-------------------------------------------------------------------------------------------------------------------
# This function is basically the Shiba transform of a given state on a bipartite lattice.
# Every singly occupied state stays unchanged. Every doubly occupied state gets transformed to a vacant state.
# Every vacant state gets transformed to a doubly occupied state.
# The entire wavefunction gets multiplied by a phase which is given by the (-1)^n, where 
# n = number of doubly occupied sites + number of singly ocupied B sites
#-------------------------------------------------------------------------------------------------------------------

def shiba(x, n):
    t = x
    i = 1
    y = 0
    count = 0
    while(i<=n):
        a = extract_lastkbits(t,2)
        if(i%2 == 0)and((a == 1)or(a == 2)):
            count = count + 1
        if(a == 3):
            count = count + 1
        if(a == 1)or(a == 2):
            y = y + (2**(2*(i - 1)))*a
        if(a == 0):
            y = y + (2**(2*(i - 1)))*3
        t = t>>2
        i = i + 1
    transf_state = (-1)**count*y
    return(transf_state)


# In[11]:


#-------------------------------------------------------------------------------------------------------------------
# This function gives out the final states for a given state after hamiltonian acts on it as an array P
# Array A gives the absolute integer values of the final states
# Array count gives the sum of coefficients of the elements in array P
# X=MS0 state, n=number of orbitals
#-------------------------------------------------------------------------------------------------------------------

import numpy as np

def final_states_transformed(x,n,t,U):
    i=1
    norm = 1/(np.sqrt(len(x)))
    P=[]
    P_coeff=[]
    #--------------------------------------
    # Creating the array P
    #--------------------------------------
    
    while(i<n):
        for j in range(0,len(x),1):
            a1=destroy_electron(x[j],i+1,'up')
            a2=create_electron(a1,i,'up')
            if(a2!=None):
                P.append(a2)
                P_coeff.append(t*norm)
            b1=destroy_electron(x[j],i+1,'down')
            b2=create_electron(b1,i,'down')
            if(b2!=None):
                P.append(b2)
                P_coeff.append(t*norm)
            c1=destroy_electron(x[j],i,'down')
            c2=create_electron(c1,i+1,'down')
            if(c2!=None):
                P.append(c2)
                P_coeff.append(t*norm)
            d1=destroy_electron(x[j],i,'up')
            d2=create_electron(d1,i+1,'up')
            if(d2!=None):
                P.append(d2)
                P_coeff.append(t*norm)
        i=i+1
            
            
    i=1
    while(i<=n):
        for j in range(0,len(x),1):
            e1=create_electron(x[j],i,'up')
            e2=create_electron(e1,i,'down')
            e3=destroy_electron(e2,i,'down')
            e4=destroy_electron(e3,i,'up')
            if(e4!=None):
                P.append(e4)
                P_coeff.append(U*norm)
        i=i+1
    
    
    b=0
    A=[]
    count=[]
    while(b<len(P)):
        j=0
        flag=0 
        #-------------------------------------------------------------------------
        # Flag variable to signal if an element has been encountered before
        # Checking if an element has occured before and updating count
        #-------------------------------------------------------------------------
        
        while(j<len(A)):
            if(abs(P[b])==A[j]):
                count[j]=count[j]+int(P[b]/abs(P[b]))*P_coeff[b]
                flag=1
            j=j+1
        
        #-------------------------------------------------------------------------
        # If the element has not been encountered, append it to the array A
        #-------------------------------------------------------------------------
        
        if(flag==0):
            A.append(abs(P[b]))
            count.append(int(P[b]/abs(P[b]))*P_coeff[b])
        b=b+1
    
#    P = np.asarray(P, dtype = int)
#    A = np.asarray(A, dtype = int)
#    count = np.asarray(count, dtype = np.float64)
        
    return(P,A,count)


# In[12]:


# Returns index of x in arr if present, else -1 
def binarySearch (arr, l, r, x):  #l and r denote the start and end indices of the array arr respectively
    # Check base case 
    if r >= l: 
  
        mid = int(l + (r - l)/2)
  
        # If element is present at the middle itself 
        if arr[mid] == x: 
            return mid 
          
        # If element is smaller than mid, then it  
        # can only be present in left subarray 
        elif arr[mid] > x: 
            return binarySearch(arr, l, mid-1, x) 
  
        # Else the element can only be present  
        # in right subarray 
        else: 
            return binarySearch(arr, mid+1, r, x) 
  
    else: 
        # Element is not present in the array 
        return -1


# In[13]:


#--------------------------------------------------------------------------------------------------------------
# This function takes input the number of orbitals as well as the number of electrons and creates 
# a compressed row format(csr) matrix storage for the Hubbard Hamiltonian
#--------------------------------------------------------------------------------------------------------------
from scipy.sparse import csr_matrix

def sparse_hamiltonian(n,e,t,U): #n=number of orbitals, e=number of electrons
    ms0states=states_with_MS0(n,e)[0]
    P_all=[]
    i=0
    while(i<len(ms0states)):
        P=final_states(ms0states[i],n,t,U)
        P_all.append(P)
        i=i+1
    
    
    row_count = np.array([0])
    col_idx_arr = [] #Stores the column index of the first unique element for the Hamiltonian
    val_arr = []
    
    L=len(P_all[0][1])
    a=0
    while(a<L):
        col_index=binarySearch(ms0states,0, len(ms0states)-1, P_all[0][1][a])
        col_val=P_all[0][2][a]
        col_idx_arr.append(col_index)
        val_arr.append(col_val)
        a=a+1
        
    
    
#    col_idx_arr = np.asarray(col_idx_arr, dtype = int)
#    val_arr  = np.asarray(val_arr, dtype = np.float64)
    
    
    
    k = 0
#    row_count.append(len(P_all[0][1])) #Stores the number of non-zero elements in a given row
    while(k < len(ms0states)):
        L=len(P_all[k][1])
        row_count = np.append(row_count, L+row_count[k])
        a=0
        if(k>0):
            while(a<L):
                col_index = binarySearch(ms0states,0, len(ms0states)-1, P_all[k][1][a])
                col_val = P_all[k][2][a]
#                col_idx_arr = np.append(col_idx_arr, col_index)
                col_idx_arr.append(col_index)
#                val_arr = np.append(val_arr, col_val)
                val_arr.append(col_val)
                a=a+1
        k=k+1
        
    dim = len(row_count) - 1
    A = csr_matrix((val_arr, col_idx_arr, row_count), shape=(dim, dim)).transpose()
    return(A)


# In[14]:


#--------------------------------------------------------------------------------------------------------------
# This function takes input the number of orbitals as well as the number of electrons and creates 
# a compressed row format(csr) matrix storage for the Hubbard Hamiltonian with the necessary symmetries
#--------------------------------------------------------------------------------------------------------------
from scipy.sparse import csr_matrix
import numpy as np

def sparse_hamiltonian_modified(n,e,t,U): #n=number of orbitals, e=number of electrons    
    ms0states=states_with_MS0(n,e)[0]
    ionic = ionic_covalent_space(n, e)[2] # The ionic space of the Hamiltonian
    covalent = ionic_covalent_space(n, e)[1] # The covalent space of the Hamiltonian
    
    first_i = [] # Stores the first element of each mixed state in the ionic space
    first_c = [] # Stores the first element of each mixed state in the covalent space
    
    for i in range(0,len(ionic),1):
        first_i.append(ionic[i][0])
    for i in range(0,len(covalent),1):
        first_c.append(covalent[i][0])
    
    row_count = np.array([0])
    col_idx_arr = [] #Stores the column index of the first unique element for the Hamiltonian
    val_arr = []
    
    
    P_all=[]
    i=0
    while(i<len(ionic)):
        P=final_states_transformed(ionic[i],n,t,U)
        P_all.append(P)
        i=i+1
    
    
    L=len(P_all[0][1])
    a=0
    while(a<L):
        scatt_states = P_all[0][1]
        scatt_coeffs = P_all[0][2]
        col_index = binarySearch(first_i,0, len(first_i)-1, scatt_states[a]) # Finds the column index of first non-zero element
        if(col_index != -1): # Proceed to next step only  if the element is found in the list first_i
            norm = np.sqrt(len(ionic[col_index])) # Finds the normalization factor of the the given element
            col_val = scatt_coeffs[a]*norm
            col_idx_arr.append(col_index)
            val_arr.append(col_val)
            shiba_el = abs(shiba(scatt_states[a], n))
            if (shiba_el != abs(scatt_states[a])):
                z = scatt_states.index(shiba_el)
                P_all[0][1].pop(z)
                P_all[0][2].pop(z)
                L = L - 1
        a=a+1
        
    
    k = 0
#    row_count.append(len(P_all[0][1])) #Stores the number of non-zero elements in a given row
    while(k < len(ionic)):
        L = len(P_all[k][1])
        a=0
        if(k>0):
            while(a<L):
                scatt_states = P_all[k][1]
                scatt_coeffs = P_all[k][2]
                col_index = binarySearch(first_i,0, len(first_i)-1, scatt_states[a])
                if(col_index != -1):
                    norm = np.sqrt(len(ionic[col_index]))
                    col_val = scatt_coeffs[a]*norm
                    col_idx_arr.append(col_index)
                    val_arr.append(col_val)
                    shiba_el = abs(shiba(scatt_states[a], n))
                    if (shiba_el != abs(scatt_states[a])):
                        z = scatt_states.index(shiba_el)
                        P_all[k][1].pop(z)
                        P_all[k][2].pop(z)
                        L = L - 1
                a=a+1
        L = len(P_all[k][1])
        row_count = np.append(row_count, len(val_arr))
        k=k+1
        
    dim = len(row_count) - 1
    A = csr_matrix((val_arr, col_idx_arr, row_count), shape=(dim, dim)).transpose()
    
    
    #----------------------------------------------------------------------------------------
    # Basically a repeat for the algorithm for the ionic subspace in the covalent subspace
    #----------------------------------------------------------------------------------------
    
    row_count = np.array([0])
    col_idx_arr = [] #Stores the column index of the first unique element for the Hamiltonian
    val_arr = []
    
    P_all=[]
    i=0
    while(i<len(covalent)):
        P=final_states_transformed(covalent[i],n,t,U)
        P_all.append(P)
        i=i+1
    
    
    L=len(P_all[0][1])
    a=0
    while(a<L):
        scatt_states = P_all[0][1]
        scatt_coeffs = P_all[0][2]
        col_index = binarySearch(first_c,0, len(first_c)-1, scatt_states[a]) # Finds the column index of first non-zero element
        if(col_index != -1): # Proceed to next step only  if the element is found in the list first_i
            norm = np.sqrt(len(covalent[col_index])) # Finds the normalization factor of the the given element
            col_val = scatt_coeffs[a]*norm
            col_idx_arr.append(col_index)
            val_arr.append(col_val)
            shiba_el = abs(shiba(scatt_states[a], n))
            if (shiba_el != abs(scatt_states[a])):
                z = scatt_states.index(shiba_el)
                P_all[0][1].pop(z)
                P_all[0][2].pop(z)
                L = L - 1
        a=a+1
        
    
    k = 0
#    row_count.append(len(P_all[0][1])) #Stores the number of non-zero elements in a given row
    while(k < len(covalent)):
        L = len(P_all[k][1])
        a=0
        if(k>0):
            while(a<L):
                scatt_states = P_all[k][1]
                scatt_coeffs = P_all[k][2]
                col_index = binarySearch(first_c,0, len(first_c)-1, scatt_states[a])
                if(col_index != -1):
                    norm = np.sqrt(len(covalent[col_index]))
                    col_val = scatt_coeffs[a]*norm
                    col_idx_arr.append(col_index)
                    val_arr.append(col_val)
                    shiba_el = abs(shiba(scatt_states[a], n))
                    if (shiba_el != abs(scatt_states[a])):
                        z = scatt_states.index(shiba_el)
                        P_all[k][1].pop(z)
                        P_all[k][2].pop(z)
                        L = L - 1
                a=a+1
        L = len(P_all[k][1])
        row_count = np.append(row_count, len(val_arr))
        k=k+1
        
    dim = len(row_count) - 1
    B = csr_matrix((val_arr, col_idx_arr, row_count), shape=(dim, dim)).transpose()
    return(A,B)

