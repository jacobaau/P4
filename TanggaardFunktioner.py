import math;
import matplotlib.pyplot as plt;
import numpy as np;
from copy import deepcopy;
import random;
#C = [0,0,0,0];
#X = [0,0,0,0];
#f = X[0]*C[0]+X[1]**2*C[1]+C[2]**X[2]+X[3]**C[3];

Variables = [];
#OChar[i] = 0 Priority, 1 Priority Add, 2 Communative, 3 IdentityLeft, 4 IdentityRight, 5 AnnihilatorLeft, 6 AnnihilatorRight, 7 Anti [0 Anti for first value, 1 Place for first value, 2 Anti for second value, 3 Place for second value],
#8 Derivatives [0 No Variables, 1 Front Variable, 2 Back Variable, 3 Both Variables]
parPr = 10;
OChar = [[0,parPr,0,None,None,None,None,[None,None,None,None],[None,None,None,None]], #0 Parenthesis start,
         [parPr,-parPr,0,None,None,None,None,[None,None,None,None],[None,None,None,None]], #1 Parenthesis end,
         [1,0,1,0,0,None,None,[3,1,3,1],[[0],[['3']],[['2']],[['2'],[2],['3']]]], #2 Add,
         [1,0,1,None,0,None,None,[2,1,2,1],[[0],[0,[2],-1,[4],['3']],[['2']],[['2'],[2],-1,[4],['3']]]], #3 Subtract,
         [2,0,1,1,1,0,0,[5,1,5,1],[[0],[['0'],[4],['3']],[['2'],[4],['1']],[['2'],[4],['1'],[2],['3'],[4],['0']]]], #4 Multiply,
         [3,0,1,None,1,0,0,[4,1,4,1],[[0],[0,[0],0,[2],-1,[4],['0'],[4],['3'],[1],0,[5],0,[0],0,[0],['1'],[1],0,[6],2,[1],0],[['2'],[5],['1']],[0,[0],['2'],[4],['1'],[2],-1,[4],['3'],[4],['0'],[1],0,[5],0,[0],0,[0],['1'],[1],0,[6],2,[1],0]]], #5 Divide,
         [4,0,0,None,1,1,None,[6,1,7,0],[[0],[['3'],[4],0,[9],['0'],[4],['0'],[6],['1']],[['2'],[4],['1'],[4],['0'],[6],0,[0],['1'],[2],-1,[4],1,[1],0],[0]]], #6 Power,
         [6,0,0,None,None,None,None,[None,None,6,0],[[0],[1,[5],['1'],[5],0,[9],['0']],[0],[0]]], #7 Log,
         [7,0,0,None,1,None,None,[None,None,9,0],[[0],[0,[0],['3'],[1],0,[4],0,[8],['1']],[0],[0]]], #8 ExponentialPower,
         [7,0,0,None,None,None,None,[None,None,8,0],[[0],[1,[5],['1']],[0],[0]]], #9 Natural Log,
         [5,0,0,None,None,None,None,[None,None,None,None],[[0],[0],[0],[0]]], #10 Factorial,
         [0,0,0,None,None,None,None,[None,None,None,None],[[['2'],[11],['3']],[['2'],[11],['3']],[['2'],[11],['3']],[['2'],[11],['3']]]] #11 Comma;
         ];
#AOChar[]

def Priority(Eq):
    if len(Eq) == 1: return[];
    Pr = [0 for i in range(0,len(Eq)//2)];
    HPar = [];
    pAdd = 0;
    i = 0;
    while i < len(Eq)//2:
        if len(HPar) > 0 and Eq[2*i+1] != [0]:
            s = OChar[Eq[2*i+1][0]][1] + pAdd - parPr;
            if OChar[Eq[2*i+1][0]][0] + s < HPar[-1] or (OChar[Eq[2*i+1][0]][0] + s <= HPar[-1] and len(Eq[2*i+1]) == 1) :
                pAdd -= parPr;
                HPar.pop(-1);
                continue;
        pAdd += OChar[Eq[2*i+1][0]][1];
        if len(Eq[2*i+1]) > 1:
            if Eq[2*i+1][1] == 0:
                HPar.append(OChar[Eq[2*i+1][0]][0] + pAdd);
                pAdd += parPr;
        Pr[i] = OChar[Eq[2*i+1][0]][0] + pAdd;
        i += 1;
    m = 0;
    if Pr[0] >= parPr:
        m = 1;
        for i in Pr[1:]:
            if i <= Pr[0]: m = 0;
    if m == 1:
        Pr = [i-parPr for i in Pr];
    return Pr;

def Operations(Eq1,Eq2,operation):
    if not isinstance(Eq1,list):
        if not isinstance(Eq2,list):
            if operation == 2:
                return Eq1 + Eq2;
            if operation == 3:
                return Eq1 - Eq2;
            if operation == 4:
                return Eq1 * Eq2;
            if operation == 5:
                return Eq1 / Eq2;
            if operation == 6:
                return Eq1 ** Eq2;
            if operation == 7:
                return math.log(Eq2, Eq1);
            if operation == 8:
                return math.exp(Eq2);
            if operation == 9:
                return math.log(Eq2);
            if operation == 10:
                return int(math.factorial(Eq1) / math.factorial(Eq2));

def RemoveParenthesis(Eq,Pr=None):
    if Pr == None:
        Pr = [0] + Priority(Eq) + [0];
    i = 1;
    while i < len(Pr)-2:
        if Eq[2*i-1] == [0]:
            j = i+1;
            low = 1000;
            add = 0;
            if 2*i-3 > 0:
                if len(Eq[2*i-3]) > 1 or Eq[2*i-3] == [5]:
                    add = 0.5;
            while Pr[i] != Pr[j]:
                low = min(low,Pr[j]-10);
                if low < Pr[i-1]+add:
                    break;
                j += 1;
            if Eq[2*j-1] == [1]:
                if low >= Pr[j+1]: 
                    Pr.pop(i);
                    Pr.pop(j-1);
                    Eq.pop(2*i-2);
                    Eq.pop(2*i-2);
                    Eq.pop(2*j-3);
                    Eq.pop(2*j-3);
                    Pr = [0] + Priority(Eq) + [0];
                    i -= 1;
        i += 1;
    return Eq;

def Evaluate(Eq):
    if Eq == None:
        return None;
    if not isinstance(Eq, list):
        return None;
    if len(Eq) == 1:
        return Eq;
    lOld = 0;
    while lOld != len(Eq) and len(Eq) > 1:
        Pr = [0] + Priority(Eq) + [0];
        Eq = RemoveParenthesis(Eq,Pr);
        Pr = [0] + Priority(Eq) + [0];
        lOld = len(Eq);
        i = 1;
        while i < len(Pr)-1:
            if len(Eq[2*i-1]) == 1:
                if Pr[i-1] <= Pr[i] and Eq[2*i-2] == OChar[Eq[2*i-1][0]][3]:
                    s = 1;
                    if 2*i-3 >= 0:
                        s = OChar[Eq[2*i-3][0]][2];
                    if Pr[i-1] < Pr[i] or s == 1:
                        Eq.pop(2*i-2);
                        Eq.pop(2*i-2);
                        Pr.pop(i);
                        i -= 1;
                elif Pr[i] >= Pr[i+1] and Eq[2*i] == OChar[Eq[2*i-1][0]][4]:
                    Eq.pop(2*i-1);
                    Eq.pop(2*i-1);
                    Pr.pop(i);
                    i -= 1;
                elif Pr[i-1] <= Pr[i] and Eq[2*i-2] == OChar[Eq[2*i-1][0]][5] and Pr[i-1] > Pr[i]:
                    print(5,Eq[2*i-1])
                    s = 1;
                    if 2*i-3 >= 0:
                        s = OChar[Eq[2*i-3][0]][2];
                    if Pr[i-1] < Pr[i] or s == 1:
                        j = i+1;
                        while Pr[i] <= Pr[j]:
                            j += 1;
                        while i != j:
                            Eq.pop(2*i-2);
                            Eq.pop(2*i-2);
                            Pr.pop(i);
                            j -= 1
                elif Pr[i] >= Pr[i+1] and Eq[2*i] == OChar[Eq[2*i-1][0]][6] and Pr[i] < Pr[i+1]:
                    j = i-1;
                    while Pr[i] <= Pr[j]:
                        j -= 1;
                    while j != i:
                        Eq.pop(2*j);
                        Eq.pop(2*j);
                        Pr.pop(j+1);
                        i -= 1;
                if i < 1: i = 1; continue;
                if Pr[i-1] <= Pr[i] and Pr[i] >= Pr[i+1]:
                    Eval = Operations(Eq[2*i-2],Eq[2*i],Eq[2*i-1][0]);
                    if Eval != None:
                        Eq[2*i-2] = Eval;
                        Eq.pop(2*i-1);
                        Eq.pop(2*i-1);
                        Pr.pop(i);
                        i -= 1;
                        i = 0;
            if len(Pr) < 3: break;
            i += 1;
    if len(Eq) == 1: return Eq;
    i = 1;
    Eq = RemoveParenthesis(Eq);
    Pr = [0] + Priority(Eq) + [0];
    while i < len(Pr)-1:
        if len(Eq[2*i-1]) > 1:
            par = 0;
            if Eq[2*i-1][1] == 0:
                j = i+1;
                while j < len(Pr)-1:
                    if Pr[j] <= Pr[i]:
                        break;
                    j += 1;
                for k in range(Eq[2*i-1][3],Eq[2*i-1][4]+1):
                    Eval = deepcopy(Eq[2*i:2*j-1]);
                    Eval = Evaluate(InsertVariables(Eval,[Eq[2*i-1][2]],[[k]]));
                    if len(Eval) > 1:
                        continue;
                    Eq[2*i-2] = Operations(Eq[2*i-2],Eval[0],Eq[2*i-1][0]);
                    if k == Eq[2*i-1][3]:
                        Eq[2*i-1][3] += 1;
                    elif k == Eq[2*i-1][4]:
                        Eq[2*i-1][4] -= 1;
                    else:
                        if par == 0:
                            par = 1;
                            Eq = Eq[:2*i-2] + [0,[0]] + Eq[2*i-2:2*j-1] + [[1],0] + Eq[2*j-1:];
                            i += 1;
                            j += 1;
                        Eq = Eq[:2*j-1] + [[Eq[2*i-1][0]]] + [0] + [Eq[2*i-1][:]] + Eq[2*i:2*j-1] + Eq[2*j-1:];
                        s = j-i;
                        Eq[2*i-1][4] = k-1;
                        if Eq[2*i-1][3] == Eq[2*i-1][4]:
                            Eval = Eq[2*i:2*j-1];
                            v = 0;
                            if len(Eval) == 1:
                                Eval = Eval[0][:];
                                v = 1;
                            Eval = InsertVariables(Eval,[Eq[2*i-1][2]],[Eq[2*i-1][3]]);
                            if v == 1: Eval = [Eval];
                            Eq = Eq[:2*i-1] + [[Eq[2*i-1][0]]] + Eval + Eq[2*j-1:];
                            i -= s-1;
                            j -= s-1;
                        i += s+1;
                        j += s+1;
                        Eq[2*i-1][3] = k+1;
                if Eq[2*i-1][3] > Eq[2*i-1][4]:
                    Eq = Eq[:2*i-1] + Eq[2*j-1:];
                    i -= 1;
        i += 1;
        Eq = RemoveParenthesis(Eq);
        Pr = [0] + Priority(Eq) + [0];
    Eq = RemoveParenthesis(Eq);
    return Eq;

def InsertVariable(Eq,Var,Value):
    i = -2;
    while i < len(Eq)-2:
        i += 2;
        if isinstance(Eq[i],list):
            if Eq[i] == Var:
                Eq = Eq[:i] + Value + Eq[i+1:];
            else:
                Eq = Eq[:i] + [InsertVariable(Eq[i],Var,Value)] + Eq[i+1:];
        if i+1 >= len(Eq): continue;
        if isinstance(Eq[i+1],list):
            if len(Eq[i+1]) > 1:
                if Eq[i+1][1] == 1:
                    if Eq[i+1][2] == Var:
                        if len(Value) == 1 and not isinstance(Value[0],list):
                            j = 0;
                            while j < len(Eq[i+1]) - 3:
                                if Eq[i+1][j+3][0] >= Value[0]:
                                    break;
                                j += 1;
                            start = i + 2;
                            end = i + 3;
                            k = i + 3;
                            w = 0;
                            while  w < len(Eq[i+1]) - 2:
                                if Eq[k] == [11]:
                                    w += 1;
                                    if j == w:
                                        start = k + 1;
                                    if j == w - 1:
                                        end = k;
                                k += 2;
                            Eq = Eq[:i+1] + [[0]] + Eq[start:end] + Eq[k:];
                    else:
                        Eq[i+1][2] = InsertVariable(Eq[i+1][2],Var,Value);
    return Eq;

def InsertVariables(Eq,Vars,Values):
    for i in range(0,len(Vars)):
        Eq = InsertVariable(Eq,Vars[i],Values[i]);
    Eq = FindVariables(Eq);
    Eq = Evaluate(Eq);
    Eq = Evaluate(Eq);
    return Eq;

def DefineVariables(Vars, Values):
    for i in range(0,len(Vars)):
        for j in range(0,len(Variables)//2):
            if Vars[i] == Variables[2*j]:
                Variables[2*j+1] == Values[i];
                Vars.pop(i);
                Values.pop(i);
                break;
    for i in range(0,len(Vars)):
        Variables.append(Vars[i]);
        Variables.append(Values[i]);

def FindVariables(Eq):
    i = -2;
    while i < len(Variables)-2:
        i += 2;
        Eq = InsertVariable(Eq,Variables[i],Variables[i+1]);
    return(Eq);
            
def Differentiate(Eq,Var):
    parenthesis = 0;
    Parsta = [0];
    while len(Eq) > 4:
        if Eq[1][0] == 0 and Eq[-2][0] == 1:
            if Parsta != [0]: break;
            i = 1;
            s = 1;
            while s > 0:
                i += 2;
                if Eq[i][0] == 0: s += 1;
                if Eq[i][0] == 1: s -= 1;
            if i == len(Eq)-2:
                if Eq[1] != [0]:
                    Parsta = Eq[1];
                parenthesis = 1;
                Eq = Eq[2:-2];
            else: break;
        else: break;
    if len(Eq) == 1:
        if Eq[0] == Var:
            return [1];
        return [0];
    Pr = Priority(Eq);
    b = 0;
    i = 0;
    while i < 10:
        j = 1;
        while j < len(Eq):
            if Pr[j//2] == i:
                if len(Eq[j]) > 1:
                    Eq = AdvancedDifferentiate(Eq,Var);
                    b = 1;
                    break;
                if OChar[Eq[j][0]][8][0] != None:
                    backVar = 0;
                    frontVar = 0;
                    k = 0;
                    while k < len(Eq[:j]):
                        if Eq[k] == Var: backVar = 1; break;
                        k += 2;
                    k = j+1;
                    while k < len(Eq):
                        if Eq[k] == Var: frontVar = 1; break;
                        k += 2;
                    Diff1 = Differentiate(Eq[:j],Var);
                    Diff2 = Differentiate(Eq[j+1:],Var);
                    if OChar[Eq[j][0]][0] > 0:
                        Diff1 = [0,[0]] + Diff1 + [[1],0];
                        Diff2 = [0,[0]] + Diff2 + [[1],0];
                    Eq = InsertVariables(OChar[Eq[j][0]][8][frontVar+backVar*2][:],[['0'],['1'],['2'],['3']],[Eq[:j],Eq[j+1:],Diff1,Diff2]);
                    b = 1;
                    break;
            j += 2;
        if b == 1: break;
        i += 1;
    if parenthesis == 1:
        return [0] + [Parsta] + Eq + [[1],0];
    return Eq;
        
def CompareVariable(Eq,Var,Change,Start,End):
    I = [1 for i in range(0,len(Change))];
    if Eq[0] == []:
        return I;
    elif len(Eq) != len(Var):
        return I;
    fail = 0;
    for i in range(0,len(Eq)):
        if fail == 1: break;
        if isinstance(Eq[i],list):
            if not isinstance(Var[i],list):
                for j in range(0,len(Change)):
                    if Eq[i] == Change[j]:
                        if not isinstance(I[j],list): I[j] = [];
                        if Start[j] <= Var[i] and Var[i] <= End[j]:
                            I[j].append(Var[i]);
                        else: fail = 1;
                        break;
                    if j == len(Change)-1: fail = 1;
            else:
                I = CompareVariable(Eq[i],Var[i],Change,Start,End);
        elif Eq[i] != Var[i]:
            fail = 1;
            break;
    if fail == 1: I = [[] for i in range(0,len(Change))];
    return I;

def AdvancedDifferentiate(Eq,Var):
    ParSta = [];
    if Eq[1][1] == 0:
        Ii = [];
        i = 2;
        Cha = [Eq[1][2]];
        Sta = [Eq[1][3]];
        End = [Eq[1][4]];
        l = len(Eq);
        while i < len(Eq):
            if isinstance(Eq[i],list):
                Ii.append(CompareVariable(Eq[i],Var,Cha,Sta,End));
            i += 2;
        I = [1 for i in range(0,len(Cha))];
        for i in range(0,len(I)):
            for j in range(0,len(Ii)):
                if not isinstance(Ii[j][i],list): break;
                if not isinstance(I[i],list): I[i] = [];
                for k in Ii[j][i]:
                    if k not in I[i]:
                        I[i].append(k);
        Il = 0;
        for i in range(0,len(I)):
            if isinstance(I[i],list):
                Il += len(I[i]);
        if Il > 0 and Il < 0:
            Len = [0 for i in Cha];
            ParSta = [0 for i in Cha];
            for i in range(0,len(I)):
                if isinstance(I[-1-i],list):
                    Len[-1-i] = len(I[-1-i]);
                else: Cha.pop(-1-i); Len.pop(-1-i);
            Pla = [1];
            for i in range(0,len(Len)-1):
                Pla.append(Pla[i]*Len[i]);
            to = Pla[-1]*Len[-1];
            Val = [0 for i in Len];
            for i in range(0,to):
                for j in range(0,len(Val)):
                    Val[j] = I[j][i//Pla[j]%Len[j]];
                Part = InsertVariables(Eq[2:],Cha,Val);
                for j in range(0,len(Cha)):
                    for k in range(0,len(Eq)//2):
                        if Eq[2*k+1][2] == Cha[j]:
                            if Val[j] == Eq[2*k+1][3]:
                                Eq[2*k+1][3] += 1;
                            elif Val[j] == Eq[2*k+1][4]:
                                Eq[2*k+1][4] -= 1;
                            else:
                                ParSta[j] = 2*k;
                                s = len(Eq[2*k:])+1;
                                Eq =  [Eq[2*k],Eq[2*k+1][:]]+Eq[2*k+2:]+[[Eq[2*k+1][0]]]+[Eq[2*k],Eq[2*k+1][:]]+Eq[2*k+2:l];
                                Eq[2*k+1][4] = Val[j]-1;
                                Eq[2*k+1+s][3] = Val[j]+1;
                            break;
                Eq = Eq + [[Eq[1][0]]] + Part;
            Eq = [0,[0]] + Differentiate(Eq,Var) + [[1],0];
        else:
            Diff = [0,[0]] + Differentiate(Eq[2:],Var) + [[1],0];
            Eq = Eq[:2] + Diff;
    for i in ParSta:
        Eq = Eq[:i] + [0,[0]] + Eq[i:] + [[1],0];
    return Eq;

def IdentifyVariable(Eq,Var):
    I = [];
    i = 0;
    while i < len(Eq):
        if Eq[i] == Var:
            I.append(i);
        i += 2;
    return I;


def IsolateVariable(Eq, AntiEq,Var,Add = []):
    if not isinstance(Eq, list):
        return None;
    if not isinstance(AntiEq, list):
        return None;
    Eq = RemoveParenthesis(Eq);
    Pr = Priority(Eq);
    I = IdentifyVariable(Eq,Var);
    if len(I) == 0:
        return [];
    b = 0;
    for i in range(1,parPr):
        j = 0;
        Itemp = I[:];
        back = 0;
        front = 1;
        while j < len(Pr):
            if 2*j in Itemp:
                back = 1;
                Itemp.remove(2*j);
                if len(Itemp) == 0:
                    front = 0;
            if Pr[j] == i:
                if len(Eq[2*j+1]) > 1:
                    AntiEq = IsolateVariable(Eq[2:], AntiEq,Var,Eq[:2]);
                    b = 1
                    break;
                if back == 0:
                    if front == 1:
                        if OChar[Eq[2*j+1][0]][7][2] == None:
                            return None;
                        if OChar[Eq[2*j+1][0]][7][3] == 1:
                            AntiEq = AntiEq + [[OChar[Eq[2*j+1][0]][7][2]]] + Add + Eq[:2*j+1];
                        else:
                            AntiEq = Eq[:2*j+1] + [[OChar[Eq[2*j+1][0]][7][2]]] + AntiEq;
                        AntiEq = [0,[0]] + AntiEq + [[1],0];
                        Eq = Eq[2*j+2:];
                        if len(Eq) == 1:
                            b = 1;
                            break;
                        AntiEq = IsolateVariable(Eq,AntiEq,Var,deepcopy(Add));
                        if AntiEq == None:
                            return None;
                        b = 1
                        break;
                if back == 1:
                    if front == 0:
                        if OChar[Eq[2*j+1][0]][7][0] == None:
                            return None;
                        if OChar[Eq[2*j+1][0]][7][1] == 1:
                            AntiEq = AntiEq + [[OChar[Eq[2*j+1][0]][7][0]]] + Add + Eq[2*j+2:];
                        else:
                            AntiEq = Eq[2*j+2:] + [[OChar[Eq[2*j+1][0]][7][0]]] + AntiEq;
                        AntiEq = [0,[0]] + AntiEq + [[1],0];
                        Eq = Eq[:2*j+1];
                        if len(Eq) == 1:
                            b = 1;
                            break;
                        AntiEq = IsolateVariable(Eq,AntiEq,Var,deepcopy(Add));
                        if AntiEq == None:
                            return None;
                        b = 1
                        break;
                    if front == 1:
                        return None;
            j += 1;
        if b == 1:
            break;
    return AntiEq;

def MinimizeVariableAnalytical(Eq,Var):
    Diff = Evaluate(Differentiate(Eq,Var));
    Step = IsolateVariable(Diff,[0],Var);
    if Step == []:
        return 1;
    Value = Evaluate(Step);
    if Value == None:
        return 0;
    DefineVariables([Var],[Value])
    return 1;

def Criteria(Points,polyDeg,differentiable):
    for i in range(0,differentiable+2):
        if len(Points) > i:
            for j in range(0,len(Points[0])):
                if j < len(Points[i]):
                    if isinstance(Points[i][j],list):
                        Points[i][j] = Points[i][j][0];
                else: Points[i].append(None);
        else:
            Points.append([]);
            for j in range(0,len(Points[0])):
                Points[i].append(None);
    A = [];
    for i in range(0,(len(Points[0])-1)*(differentiable+1)*2):
        vPoly = i//((differentiable+1)*2);
        diff = (i%((differentiable+1)*2))//2;
        end = i%2;
        if((vPoly+end == len(Points[0])-1) and Points[1+diff][vPoly+end] == None): continue;
        if(end == 0 and Points[1+diff][vPoly] == None): continue;
        A.append([]);
        for j in range(0,(len(Points[0])-1)*(polyDeg+1)+1):
            hPoly = j//(polyDeg+1);
            deg = j%(polyDeg+1);
            if hPoly == len(Points[0])-1:
                if Points[1+diff][vPoly+end] == None: A[len(A)-1].append(0);
                else: A[len(A)-1].append(Points[1+diff][vPoly+end]);
            else:
                if (vPoly == hPoly):
                    if diff <= deg:
                        A[len(A)-1].append((math.factorial(deg)/math.factorial(deg-diff))*Points[0][vPoly+end]**(deg-diff));
                    else: A[len(A)-1].append(0);
                elif (vPoly+end == hPoly and Points[1+diff][vPoly+end] == None):
                    if diff <= deg:
                        A[len(A)-1].append(-(math.factorial(deg)/math.factorial(deg-diff))*Points[0][vPoly+end]**(deg-diff));
                    else: A[len(A)-1].append(0);
                else: A[len(A)-1].append(0);
    return A;


#Gauss Elimination
def GaussEli(B):
    Pivots= [None for i in range(0,len(B[0])-1)];
    for i in range(0,len(B[0])-1):
        for j in range(0,len(B)):
            if B[j][i]!=0:
                Pivots[i]=j;
                for k in range(0,i):
                    if B[j][k]!=0: Pivots[i]=None;
            if Pivots[i]!=None:
                for k in range(0,len(B)):
                    if Pivots[i]==k:
                        div = B[Pivots[i]][i];
                        for l in range(0,len(B[0])):
                            B[k][l] = B[k][l]/div;
                    else:
                        
                        mul = B[k][i]/B[Pivots[i]][i];
                        for l in range(0,len(B[0])):
                            B[k][l] = B[k][l]-B[Pivots[i]][l]*mul;
                break;
    return [B,Pivots];

#Substitute coefficients with variables
def SubVar(Matrix,Pivots,deg):
    n = len(Pivots);
    Variables = [];
    for i in range(0,n):
        if Pivots[i] == None: Variables.append(i);
    Subs = [];
    for i in range(0,len(Matrix[0])-1):
        Subs.append([]);
        if Pivots[i] == None: Subs[i].append([i//(deg+1),[11],i%(deg+1)]); continue;
        Subs[i].append(Matrix[Pivots[i]][-1]);
        for j in Variables:
            if Matrix[Pivots[i]][j] == 0: continue;
            Subs[i].append([2]);
            Subs[i].append(-Matrix[Pivots[i]][j]);
            Subs[i].append([4]);
            Subs[i].append([j//(deg+1),[11],j%(deg+1)]);
    return Subs;

def GenerateRandom(lists,amount,minimum,maximum):
    R = [];
    for i in range(0,lists*amount):
        R.append(random.uniform(minimum,maximum));
    return R;

def GenerateVariableNames(First,amount):
    Vars = [];
    for i in First:
        for j in range(1,amount+1):
            Vars.append([i,[11],j]);
    return Vars;

def NelsonSiegelForm(terms,Variable=[-1,[11],[1]]):
    Eq = [[0,[11],0],[2],[0,[11],1],[4],0,[0],1,[2],-1,[4],0,[8],0,[0],0,[2],-1,[4]]+[Variable]+[[5],[0,[11],-1],[1],0,[1],0,[5]]+[Variable]+[[4],[0,[11],-1]]
    Eq = Eq + [[2],[0,[11],2],[4],0,[0],0,[0],1,[2],-1,[4],0,[8],0,[0],0,[2],-1,[4]]+[Variable]+[[5],[0,[11],-1],[1],0,[1],0,[5]]+[Variable]+[[4],[0,[11],-1],[2],-1,[4],0,[8],0,[0],0,[2],-1,[4]]+[Variable]+[[5],[0,[11],-1],[1],0,[1],0]
    return(Eq);

def Polynomial(deg,var,X=[-1,[11],[1]]):
    Eq = [[var,[11],0]]; 
    for i in range(1,deg+1):
        if i == 1:
            Eq = Eq + [[2],[var,[11],i],[4],X];
        else:
            Eq = Eq + [[2],[var,[11],i],[4],X,[6],i];
    return Eq;

def PieceWisePoly(deg,Coefs,Def,X=[-1,[11],[1]]):
    Eq = [0,[0,1,X] + Def]
    M,P = GaussEli(Criteria([[[0]]+Def+[[30]]],deg,2))
    Subs = SubVar(M,P,deg);
    for i in Coefs:
        Eq = Eq + Polynomial(deg,i[0],X) + [[11]];
    Eq = Eq + [0,[1],0];
    Vars = [[Coefs[i//(deg+1)][0],[11],i%(deg+1)] for i in range(0,(deg+1)*len(Coefs))]
    i = 0;
    DefSubs = [];
    DefVars = [];
    while i < len(Subs):
        if Subs[i][0] != Vars[i]:
            DefSubs.append([0,[0]] + Subs[i] + [[1],0]);
            DefVars.append(Vars[i]);
            Subs.pop(i);
            Vars.pop(i);
            i -= 1;
        i += 1;
    DefineVariables(DefVars,DefSubs);
    Eq = FindVariables(Eq);
    return Eq, Vars;

def RSS(Eq,start,end,Y=[-2,[11],[1]]):
    return [0,[2,0,[1],start,end],0,[0]]+Eq+[[2],-1,[4]]+[Y]+[[1],0,[6],2];

def Smooth(Eq,wei,Mat,Vars):
    Eq.append([2]);
    Eq.append(wei);
    Eq.append([4]);
    Eq.append(0);
    Eq.append([0]);
    Eq.append(0);
    for i in range(0,len(Mat)):
        Eq.append([2]);
        Eq.append(Mat[i]);
        Eq.append([4]);
        Eq.append(Vars[i//len(Vars)]);
        Eq.append([4]);
        Eq.append(Vars[i%len(Vars)]);
    Eq.append([1]);
    Eq.append(0);
    return Eq;




