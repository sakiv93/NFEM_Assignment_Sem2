''' Author: Bauskar Nikhil '''
''' eMail: nikhil.bauskar@student.tu-freiberg.de '''
import numpy as np

def materialRoutine(epsilon,epsilon_pl, T_m):
    lamda = 14.0
    mu = 28.0
    sigma_y = 70.0
    hardening_e_pl = 0.0

    C_el = np.array([[lamda+2*mu, lamda, lamda],
                    [lamda, lamda+2*mu, lamda],
                    [lamda, lamda, lamda+2*mu]
                    ]) 
    
    sigma_trial = np.matmul(C_el,(epsilon-epsilon_pl))             # trial stress
    sigma_dev_trial = sigma_trial - 1/3 * (np.sum(sigma_trial))
    sigma_equivalent_trial = (1.5*np.sum(sigma_dev_trial*sigma_dev_trial))**0.5
    
    print(sigma_equivalent_trial,sigma_y)

    #delta_lamda = sigma_equivalent_trial/(3*mu + sigma_y)
    if sigma_equivalent_trial > sigma_y:
        delta_lamda = sigma_equivalent_trial/(3*mu + sigma_y)
    elif sigma_equivalent_trial <= sigma_y:
        delta_lamda = 0.

    #delta_lamda = sigma_equivalent_trial/(3*mu + sigma_y)
    if sigma_equivalent_trial==0:
        sigma_updated=np.zeros([3,1])
        Ct_1 = (3*lamda+2*mu) / 3
        Ct_2=0
        Ct_2_nd=0
        Ct_3=0
    else:
        if sigma_equivalent_trial > sigma_y:
            delta_lamda = sigma_equivalent_trial/(3*mu + sigma_y)
        elif sigma_equivalent_trial <= sigma_y:
            delta_lamda = 0.

    

        sigma_updated = (1/3 * np.ones([3,1]) * np.sum(sigma_trial)) + sigma_dev_trial*(1-(3*mu*delta_lamda/sigma_equivalent_trial))

        Ct_1 = (3*lamda+2*mu) / 3
        Ct_2 = (4*mu/3) * (1-(3*mu*delta_lamda/sigma_equivalent_trial))
        Ct_2_nd = (-1*mu/3) * (1-(3*mu*delta_lamda/sigma_equivalent_trial))
        Ct_3 = ((3*mu) / (sigma_equivalent_trial**2))* np.matmul(sigma_dev_trial,np.transpose(sigma_dev_trial))



    Ct = np.array([[Ct_1 + Ct_2, Ct_2_nd, Ct_2_nd],
                    [Ct_2_nd, Ct_1 + Ct_2, Ct_2_nd],
                    [Ct_2_nd, Ct_2_nd, Ct_1 + Ct_2]
                    ]) 
    Ct_out = Ct + Ct_3 

    epsilon_pl = epsilon - np.transpose(np.matmul(np.transpose(sigma_updated), np.linalg.inv(C_el)))

    return Ct_out, sigma_updated, epsilon_pl        
    
#     with np.errstate(divide='ignore'):
#         sigma_updated = (1/3 * np.ones([3,1]) * np.sum(sigma_trial)) + sigma_dev_trial*(1-(3*mu*delta_lamda/sigma_equivalent_trial))

#     Ct_1 = (3*lamda+2*mu) / 3
#     Ct_2 = (4*mu/3) * (1-(3*mu*delta_lamda/sigma_equivalent_trial))
#     Ct_2_nd = (-1*mu/3) * (1-(3*mu*delta_lamda/sigma_equivalent_trial))
#     Ct_3 = ((3*mu) / (sigma_equivalent_trial**2))* np.matmul(sigma_dev_trial,np.transpose(sigma_dev_trial))



#     Ct = np.array([[Ct_1 + Ct_2, Ct_2_nd, Ct_2_nd],
#                     [Ct_2_nd, Ct_1 + Ct_2, Ct_2_nd],
#                     [Ct_2_nd, Ct_2_nd, Ct_1 + Ct_2]
#                     ]) 
#     Ct_out = Ct + Ct_3 

#     epsilon_pl = epsilon - np.transpose(np.matmul(np.transpose(sigma_updated), np.linalg.inv(C_el)))

#     return Ct_out, sigma_updated, epsilon_pl

# #print(materialRoutine(np.array([[1],[0],[2]]),np.array([[0],[0],[0]]), 0.1))