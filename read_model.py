def fill(array, list):
    for i in list:
        array.append(float(i))

def readConstrains(prices, costs, plan, restrictions, p):
    f = open("C:/Users/Markk/Downloads/Telegram Desktop/Simplex/requirements.txt")
    for i in f.readlines():
        temp = i.replace(' ', '').split(':')
        if i[0] == '#' or len(i) == 1:
            continue
        if temp[0] == 'Prices':
            fill(prices, temp[1].split(','))
        if temp[0] == 'Costs':
            fill(costs, temp[1].split(','))
        if temp[0] == 'Plan':
            fill(plan, temp[1].split(','))
        if temp[0] == 'Restrictions':
            fill(restrictions, temp[1].split(','))
        if temp[0] == 'Product':
            p.append([])
            fill(p[-1], temp[1].split(','))

def fillMatrices():
    prices = []
    plan = []
    restrictions = []
    p = [] 
    costs = []
    
    readConstrains(prices, costs, plan, restrictions, p)
    if len(restrictions) != len(costs) or len(prices) != len(p) or len(prices) != len(plan) or any(len(i) != len(restrictions) for i in p):
        raise Exception("Model is incorrect")
    
    lessOrEqual = len(restrictions)
    moreOrEqual = len(plan)
    prices = [prices[i] - costs[0] * p[i][0] - costs[1] * p[i][1] for i in range(len(prices))]
    A = []
    c = [-i for i in prices]
    b = restrictions + plan
    for i in range(lessOrEqual):
        A.append([p[j][i] for j in range(len(p))])
    for i in range(moreOrEqual):
        A.append([0 for _ in range(len(p))])
        A[i + lessOrEqual][i] = 1
    return A, b, c, moreOrEqual, lessOrEqual
