import matplotlib.pyplot as plt
import math

#Syarat Batas
bataserror = 1e-12
Npartisi = int(input('Masukkan Jumlah Partisi\t: '));
rmax = float(input('Masukkan r maksimum\t:'));
#V0 = float(input('Masukkan Besar Potensial Barrier\t:'));
V0= 0;
BKi = float(input('Masukkan Batas Kiri Potensial Barrier\t:'));
#BKi = 0;
BKa = float(input('Masukkan Batas Kanan Potensial Barrier\t:'));
#BKa = 0;
dx = rmax/Npartisi;
l=0;
h=1;
m=1/2;
bataskiri = 0;
bataskanan = 0;

#Potential energy function
def V(r):
        #t=math.exp((r-7)/0.6);
        #return l*(l+1)/(r*r)-50*(1-5*t/(3*(1+t)))/(1+t);
        if r>=BKi and r<=BKa :
                return V0;
        return 0;

def f(psi,r,E):
        #return psi*(V(r)-E-1/(r*r)*l*(l+1));
        return psi*(V(r)-E)
def ff(psi,r,E):
        #return psi*(V(r)-E-1/(r*r)*l*(l+1));
        return psi*(-E)

#inisialisasi program
r = [(i+1)*dx for i in range (Npartisi+1)]
Psi = [0 for i in range (Npartisi+1)]
PsiSampel = [0 for i in range (Npartisi+1)]
Psi2= [0 for i in range (Npartisi+1)]
Pot = [0 for i in range (Npartisi+1)]

for i in range (Npartisi+1): #inisialisai potensial
    Pot[i]=V(r[i]);

Psi[0] = bataskiri;
PsiSampel[0] = bataskiri;
Psi2[0] = dx;

dE = 1;

#Differential Equation Solver, returning end value
#Heun Methods
def PDB(Psi,Psi2,E):
    for i in range (1,Npartisi+1):

        k1= dx * f(Psi[i-1],r[i-1],E);
        l1= dx * Psi2[i-1];
        
        k2= dx * f(Psi[i-1]+l1/2 , r[i-1]+ dx/2, E);
        l2= dx * (Psi2[i-1]+k1/2);
        
        
        k3= dx * f(Psi[i-1]+l2/2 , r[i-1]+ dx/2, E);
        l3= dx * (Psi2[i-1]+k2/2);
        
        k4= dx * f(Psi[i-1]+l3 , r[i-1]+ dx, E);
        l4= dx * (Psi2[i-1]+k3);
        

        Psi2[i]=Psi2[i-1]+1/6*(k1+2*k2+2*k3+k4);
        Psi[i]= Psi[i-1] +1/6*(l1+2*l2+2*l3+l4);
    return Psi

def PDBSampel(Psi,Psi2,E):
    for i in range (1,Npartisi+1):

        k1= dx * ff(Psi[i-1],r[i-1],E);
        l1= dx * Psi2[i-1];
        
        k2= dx * ff(Psi[i-1]+l1/2 , r[i-1]+ dx/2, E);
        l2= dx * (Psi2[i-1]+k1/2);
        
        
        k3= dx * ff(Psi[i-1]+l2/2 , r[i-1]+ dx/2, E);
        l3= dx * (Psi2[i-1]+k2/2);
        
        k4= dx * ff(Psi[i-1]+l3 , r[i-1]+ dx, E);
        l4= dx * (Psi2[i-1]+k3);
        

        Psi2[i]=Psi2[i-1]+1/6*(k1+2*k2+2*k3+k4);
        Psi[i]= Psi[i-1] +1/6*(l1+2*l2+2*l3+l4);
    return Psi       

                
def CariSimpul (Psi):
        for i in range (1,Npartisi+1):
                if r[i]>BKa:
                        if Psi[i]<1e-3 and Psi[i]>-1e-3:
                                return r[i];

def CariJumlahSimpul (Psi):
        NSimpul=0;
        for i in range (1,Npartisi+1):
                if Psi[i]*Psi[i-1]<0:
                        NSimpul=NSimpul+1;
                        if r[i]>BKa:
                                return NSimpul;
                        
def CariSimpulkeN (Psi,NCari):
        NSimpul=0;
        for i in range (1,Npartisi+1):
                if Psi[i]*Psi[i-1]<0:
                        NSimpul=NSimpul+1;
                        if NSimpul==NCari:
                                return r[i];

def Lebarantarsimpul (Psi):
        for i in range (1,Npartisi+1):
                        if Psi[i]<1e-3 and Psi[i]>-1e-3:
                                return r[i];


#Utama
plt.figure();
lgndlbl=[];
for Epot in range (10,60,10):
        Eplot = [];
        BedaFasa = [];
        V0=Epot;
        for Ei in range (1,101):
                E=Ei;
                PsiSampel = PDBSampel(PsiSampel,Psi2,E);
                Psi = PDB(Psi,Psi2,E);
                NCari = CariJumlahSimpul (PsiSampel)
                SimpulSampel = CariSimpulkeN (PsiSampel,NCari);
                Simpul = CariSimpulkeN (Psi,NCari);
                #print(E,Simpul,SimpulSampel)
                Bedafasa = (Simpul-SimpulSampel)*math.sqrt(E);
                Eplot.append (E);
                BedaFasa.append(Bedafasa);
        lgndlbl.append("E="+str(V0));
        plt.plot(Eplot,BedaFasa);
plt.title("Variasi Energi");
plt.ylabel('Phase Shift (rad)')
plt.xlabel('Energy')
plt.legend(lgndlbl);
plt.show();        

#Psi = PDB(Psi,Psi2,0);
#plt.plot(r,Psi);
#plt.show();
#plt.figure();
#plt.title("Simpul Sampel");
#plt.plot(r,PsiSampel);
#plt.figure();
#plt.title("Simpul");
#plt.plot(r,Psi);
#plt.show();


