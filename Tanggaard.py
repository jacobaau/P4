from TanggaardFunktioner import *;

#Denne funktion er ikke vigtig at forstå, længere nede er dataene og splinet

def MinimizeVariables(Eq,Vars,tolerance,iterations,smoothing,Bonds,typ):
    Numerical = Vars;
    Estimate = [[0] for i in Numerical];
    Diff1 = [];
    Smo = Smooth([0],smoothing,SmoMat,Coefficients);
    SDiff1 = [];
    for i in Numerical:
        A = deepcopy(Eq);
        S = deepcopy(Smo);
        Diff1.append(Differentiate(A,i));
        SDiff1.append(Differentiate(S,i));
    Diff2 = [];
    SDiff2 = [];
    for i in range(0,len(Numerical)):
        for j in Numerical:
            A = deepcopy(Diff1[i]);
            S = deepcopy(SDiff1[i]);
            SDiff2.append(Differentiate(S,j));
    error = tolerance + 1;
    n = 0;
    while error > tolerance and n < iterations:
        n += 1;
        D = [];
        for i in range(0,len(Numerical)):
            D.append([]);
            s = 0;
            e = 0;
            di = 0;
            dj = 0;
            for j in range(0,len(Numerical)):
                s = 0;
                for k in range(0,len(Bonds)):
                    be = 0;
                    bi = 0;
                    bj = 0;
                    bij = 0;
                    for ly in range(1,Bonds[k][3]-2018):
                        l = ly - Bonds[i][4];
                        e = InsertVariables(deepcopy(Eq),[[-1,[11],[1]]]+Numerical,[[l]]+Estimate)[0];
                        di = InsertVariables(deepcopy(Diff1[i]),[[-1,[11],[1]]]+Numerical,[[l]]+Estimate)[0];
                        dj = InsertVariables(deepcopy(Diff1[j]),[[-1,[11],[1]]]+Numerical,[[l]]+Estimate)[0];
                        dij = 0
                        if typ == 0:
                            bi += Bonds[k][2]*di;
                            bj += Bonds[k][2]*dj;
                        if typ == 1:
                            be = Bonds[k][2]*math.exp(-l*e);
                            bi += Bonds[k][2]*di*(-l)*math.exp(-l*e);
                            bj += Bonds[k][2]*dj*(-l)*math.exp(-l*e);
                            bij += Bonds[k][2]*dij*(-l)*math.exp(-l*e)+Bonds[k][2]*di*dj*l*l*math.exp(-l*e);
                        if typ == 2:
                            be = Bonds[k][2]*math.exp(-l*e/(l+1));
                            bi += Bonds[k][2]*di*(-l/(l+1))*math.exp(-l*e/(l+1));
                            bj += Bonds[k][2]*dj*(-l/(l+1))*math.exp(-l*e/(l+1));
                            bij += Bonds[k][2]*dij*(-l/(l+1))*math.exp(-l*e/(l+1))+Bonds[k][2]*di*dj*l/(l+1)*l/(l+1)*math.exp(-l*e/(l+1));
                    li = Bonds[k][3]-2018 - Bonds[i][4];
                    if typ == 0:
                        bi += Bonds[k][0]*di;
                        bj += Bonds[k][0]*dj;
                    if typ == 1:
                        be += Bonds[k][0]*math.exp(-li*e);
                        bi += Bonds[k][0]*di*(-li)*math.exp(-li*e);
                        bj += Bonds[k][0]*dj*(-li)*math.exp(-li*e);
                        bij += Bonds[k][0]*dij*(-li)*math.exp(-li*e)+Bonds[k][0]*di*dj*li*li*math.exp(-li*e);
                    if typ == 2:
                        be += Bonds[k][0]*math.exp(-li*e/(li+1));
                        bi += Bonds[k][0]*di*(-li/(li+1))*math.exp(-li*e/(li+1));
                        bj += Bonds[k][0]*dj*(-li/(li+1))*math.exp(-li*e/(li+1));
                        bij += Bonds[k][0]*dij*(-li/(li+1))*math.exp(-li*e/(li+1))+Bonds[k][0]*di*dj*li/(li+1)*li/(li+1)*math.exp(-li*e/(li+1));
                    s += bi*bj
                s = s/len(Bonds) + InsertVariables(deepcopy(SDiff2[i*len(Numerical)+j]),Numerical,Estimate)[0];
                #print(1,s);
                D[-1].append(s);
            s = 0;
            for k in range(0,len(Bonds)):
                be = 0;
                bi = 0;
                for ly in range(1,Bonds[k][3]-2018):
                    l = ly - Bonds[i][4];
                    e = InsertVariables(deepcopy(Eq),[[-1,[11],[1]]]+Numerical,[[l]]+Estimate)[0];
                    di = InsertVariables(deepcopy(Diff1[i]),[[-1,[11],[1]]]+Numerical,[[l]]+Estimate)[0];
                    if typ == 0:
                        be += Bonds[k][2]*e;
                        bi += Bonds[k][2]*di;
                    if typ == 1:
                        be += Bonds[k][2]*math.exp(-l*e);
                        bi += Bonds[k][2]*di*(-l)*math.exp(-l*e);
                    if typ == 2:
                        be += Bonds[k][2]*math.exp(-l*e/(l+1));
                        bi += Bonds[k][2]*di*(-l/(l+1))*math.exp(-l*e/(l+1));
                li = Bonds[k][3]-2018 - Bonds[i][4];
                if typ == 0:
                    be += Bonds[k][0]*e;
                    bi += Bonds[k][0]*di;
                if typ == 1:
                    be += Bonds[k][0]*math.exp(-li*e);
                    bi += Bonds[k][0]*di*(-li)*math.exp(-li*e);
                if typ == 2:
                    be += Bonds[k][0]*math.exp(-li*e/(li+1));
                    bi += Bonds[k][0]*di*(-li/(li+1))*math.exp(-li*e/(li+1));
                s += bi*(be-Bonds[k][1]);
            s = s/len(Bonds) + InsertVariables(deepcopy(SDiff1[i]),Numerical,Estimate)[0];
            #print(2,s);
            D[-1].append(s);
        D = GaussEli(D)[0];
        for i in range(0,len(D)):
            Estimate[i][0] = Estimate[i][0] - D[i][-1];
        #print(n,Estimate);
    print(n,Estimate);
    return Estimate;

#Her er dataene og splines defineret

                
#udkommenter 'Bonds' for at bruge datasættet
#2020
Bonds = [[100, 106.925, 1.5, 2023, 0.39], [100, 105.485, 3.0, 2021, 0.39], [100, 134.55, 7.0, 2024, 0.39], [100, 102.04, 0.25, 2022, 0.39], [100, 106.811, 0.5, 2029, 0.39], [100, 103.7, 0.1, 2023, 0.39], [100, 111.45, 0.1, 2030, 0.39], [100, 100.426, 0.25, 2020, 0.39], [100, 112.276, 1.75, 2025, 0.39], [100, 106.15, 0.5, 2027, 0.39], [100, 102.104, 0.25, 2052, 0.39], [100, 102.45, 0.63, 2023, 0.81], [100, 100.696, 1.0, 2020, 0.81], [100, 101.342, 0.5, 2021, 0.81], [100, 187.4, 4.5, 2039, 0.39]];

#2019
#Bonds = [[100, 109.169, 1.5, 2023, 0.39], [100, 109.114, 3.0, 2021, 0.39], [100, 141.57, 7.0, 2024, 0.39], [100, 102.92, 0.25, 2022, 0.39], [100, 105.3, 0.5, 2029, 0.39], [100, 107.48, 0.1, 2023, 0.39], [100, 113.56, 0.1, 2030, 0.39], [100, 101.356, 0.25, 2020, 0.39], [100, 113.888, 1.75, 2025, 0.39], [100, 105.74, 0.5, 2027, 0.39], [100, 102.45, 0.63, 2023, 0.81], [100, 102.02, 1.0, 2020, 0.81], [100, 101.56, 0.5, 2021, 0.81], [100, 181.53, 4.5, 2039, 0.39]];

for i in range(0,len(Bonds)):
    Bonds[i][2]=Bonds[i][2]/100*Bonds[i][0];
    #Bonds[i][3]+=1; #Udkommenter for 2019 Bonds
    Bonds[i][1]+=Bonds[i][2]*Bonds[i][4];

#Glathedsmatricen
SmoMat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4/15, 1/20, 0, 0, 0, 0, 1/20, 1/6, 5/132, 0, 0, 0, 0, 5/132, 4/33, 4/143, 0, 0, 0, 0, 4/143, 4/39]

#Spline
Equation = [0,[0,1,[-1,[11],[1]],[1],[3],[6],[11],[17],[24]],[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],1/30,[4],[-1,[11],[1]],[6],3,[2],-1/10,[4],[-1,[11],[1]],[6],2,[2],1/10,[4],[-1,[11],[1]],[2],-1/30,[1],0,[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],-1/45,[4],[-1,[11],[1]],[6],3,[2],2/5,[4],[-1,[11],[1]],[6],2,[2],-7/5,[4],[-1,[11],[1]],[2],22/15,[1],0,[2],[1,[11],3],[4],0,[0],1/72,[4],[-1,[11],[1]],[6],3,[2],-1/8,[4],[-1,[11],[1]],[6],2,[2],3/8,[4],[-1,[11],[1]],[2],-3/8,[1],0,[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],[-1,[11],[1]],[2],-10/3,[1],0,[2],[1,[11],3],[4],0,[0],-1/120,[4],[-1,[11],[1]],[6],3,[2],11/40,[4],[-1,[11],[1]],[6],2,[2],-81/40,[4],[-1,[11],[1]],[2],177/40,[1],0,[2],[1,[11],4],[4],0,[0],1/165,[4],[-1,[11],[1]],[6],3,[2],-6/55,[4],[-1,[11],[1]],[6],2,[2],36/55,[4],[-1,[11],[1]],[2],-72/55,[1],0,[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],[-1,[11],[1]],[2],-10/3,[1],0,[2],[1,[11],3],[4],0,[0],[-1,[11],[1]],[2],-20/3,[1],0,[2],[1,[11],4],[4],0,[0],-1/198,[4],[-1,[11],[1]],[6],3,[2],17/66,[4],[-1,[11],[1]],[6],2,[2],-223/66,[4],[-1,[11],[1]],[2],2669/198,[1],0,[2],[1,[11],5],[4],0,[0],1/234,[4],[-1,[11],[1]],[6],3,[2],-11/78,[4],[-1,[11],[1]],[6],2,[2],121/78,[4],[-1,[11],[1]],[2],-1331/234,[1],0,[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],[-1,[11],[1]],[2],-10/3,[1],0,[2],[1,[11],3],[4],0,[0],[-1,[11],[1]],[2],-20/3,[1],0,[2],[1,[11],4],[4],0,[0],[-1,[11],[1]],[2],-34/3,[1],0,[2],[1,[11],5],[4],0,[0],-1/273,[4],[-1,[11],[1]],[6],3,[2],24/91,[4],[-1,[11],[1]],[6],2,[2],-485/91,[4],[-1,[11],[1]],[2],9092/273,[1],0,[11],
[1,[11],0],[2],[1,[11],1],[4],[-1,[11],[1]],[2],[1,[11],2],[4],0,[0],[-1,[11],[1]],[2],-10/3,[1],0,[2],[1,[11],3],[4],0,[0],[-1,[11],[1]],[2],-20/3,[1],0,[2],[1,[11],4],[4],0,[0],[-1,[11],[1]],[2],-34/3,[1],0,[2],[1,[11],5],[4],0,[0],[-1,[11],[1]],[2],-52/3,[1],0,[11],0,[1],0];
#Koefficienter c_i
Coefficients = [[1,[11],0],[1,[11],1],[1,[11],2],[1,[11],3],[1,[11],4],[1,[11],5]];


NS = deepcopy(Equation)

Sol = [];
#MinmizeVariables(Spline,Koefficienter,ignorer,iterationer,\lambda,Data,forslag)
Sol.append(MinimizeVariables(Equation,Coefficients,0.1,1,1000,Bonds,0));#d(t)=\phi(t)
Sol.append(MinimizeVariables(Equation,Coefficients,0.1,3,1000,Bonds,1));#d(t)=e^{-t*\phi(t)}
Sol.append(MinimizeVariables(Equation,Coefficients,0.1,3,1000,Bonds,2));#d(t)=e^{-t*\phi(t)/(t+1)}

#Indsætter data til kurve
X = [i/5 for i in range(1,166)];
Y = [];
for j in range(0,len(Sol)):
    Y.append([]);
    for i in X:
        s = InsertVariables(deepcopy(NS),Coefficients+[[-1,[11],[1]]],Sol[j]+[[i]])[0];
        if j == 0:
            if s > 0:
                s = -math.log(s)/i;
            else: s = 0;
            Y[-1].append(s);
        if j == 1:
            s = s;
            Y[-1].append(s);
        if j == 2:
            s = s/(i+1);
            Y[-1].append(s);
    
plt.plot(X,Y[0],'r',label=r'$d(t) = \phi(t)$');
plt.plot(X,Y[1],'g--',label=r'$d(t) = exp(-t\cdot\phi(t))$');
plt.plot(X,Y[2],'b-.',label=r'$d(t) = exp(-t\cdot\frac{\phi(t)}{t+1})$');
plt.axvline(1,ls=':',c='k',lw = 0.5, label='knudepunkter');
for i in [3,6,11,17,24]:
    plt.axvline(i,ls=':',c='k',lw = 0.5);
Xp = [1.054,1.484,2.084,2.482,3.483,4.096,4.482,5.470,6.483,8.482,10.483,11.483,20.482];
Yp = [-0.00892,-0.00655,-0.00242,-0.00622,-0.00577,0.00036,-0.01013,-0.00418,-0.00348,-0.00165,0.00004,-0.01011,0.00845];
plt.scatter(Xp,Yp,edgecolor='k',color='none',label='Observerede nulkuponrenter');
plt.legend(loc="lower right");
plt.title(r'23 maj 2019,    $\lambda = 1000$');
plt.xlabel('Tid til udløbsdato (år)');
plt.ylabel('Nulkuponrente');
plt.show();






